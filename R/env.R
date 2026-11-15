#' Analyze environmental factors across groups
#'
#' Environmental values are taken from `data$env` and matched to
#' `data$metadata` by sample row names. By default all numeric environmental
#' variables are analyzed.
#'
#' @param data An `aq_data` object containing `env` and `metadata`.
#' @param group One metadata column defining the groups.
#' @param variables Optional character vector of environmental variables.
#'   `NULL` selects all numeric columns in `data$env`.
#' @param test Overall test: `kruskal`, `anova`, or `none`.
#' @param posthoc Pairwise test: `dunn`, `wilcox`, `both`, or `none`.
#' @param p_adjust P-value adjustment method, e.g. `BH`, `holm`, or `bonferroni`.
#' @return An `aq_env` object containing long-format data, overall tests,
#'   pairwise tests, and analysis settings.
#' @export
aq_env <- function(data, group, variables = NULL,
                   test = c("kruskal", "anova", "none"),
                   posthoc = c("both", "dunn", "wilcox", "none"),
                   p_adjust = "BH") {
  if (!inherits(data, "aq_data")) {
    stop("`data` must be an aq_data object.", call. = FALSE)
  }
  if (is.null(data$env)) {
    stop("No environmental data found. Supply env when creating aq_data().",
         call. = FALSE)
  }
  if (is.null(data$metadata)) {
    stop("Metadata is required to obtain the grouping variable.",
         call. = FALSE)
  }
  if (missing(group) || !is.character(group) || length(group) != 1L ||
      is.na(group) || !nzchar(group)) {
    stop("Please supply one metadata column through `group`.",
         call. = FALSE)
  }
  if (!group %in% names(data$metadata)) {
    stop("Grouping variable `", group, "` was not found in data$metadata.",
         call. = FALSE)
  }

  test <- match.arg(test)
  posthoc <- match.arg(posthoc)

  valid_adjust <- c("holm", "hochberg", "hommel", "bonferroni",
                    "BH", "BY", "fdr", "none")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% valid_adjust) {
    stop("`p_adjust` must be one of: ",
         paste(valid_adjust, collapse = ", "), call. = FALSE)
  }

  env <- as.data.frame(data$env, check.names = FALSE)
  if (is.null(rownames(env)) ||
      !all(rownames(data$otu) %in% rownames(env))) {
    stop("Environmental data row names must contain every OTU sample ID.",
         call. = FALSE)
  }
  env <- env[rownames(data$otu), , drop = FALSE]

  if (is.null(variables)) {
    variables <- names(env)[vapply(env, is.numeric, logical(1))]
  } else {
    if (!is.character(variables) || !length(variables) ||
        anyNA(variables) || any(!nzchar(variables))) {
      stop("`variables` must be NULL or a non-empty character vector.",
           call. = FALSE)
    }
    missing_vars <- setdiff(variables, names(env))
    if (length(missing_vars)) {
      stop("Unknown environmental variable(s): ",
           paste(missing_vars, collapse = ", "), call. = FALSE)
    }
    non_numeric <- variables[!vapply(env[variables], is.numeric, logical(1))]
    if (length(non_numeric)) {
      stop("Selected environmental variable(s) must be numeric: ",
           paste(non_numeric, collapse = ", "), call. = FALSE)
    }
  }
  if (!length(variables)) {
    stop("No numeric environmental variables were found in data$env.",
         call. = FALSE)
  }
  variables <- unique(variables)

  md <- data$metadata[rownames(env), , drop = FALSE]
  grp <- md[[group]]
  if (anyNA(grp)) {
    stop("Grouping variable `", group, "` contains missing values.",
         call. = FALSE)
  }
  grp <- droplevels(factor(grp))
  if (nlevels(grp) < 2L) {
    stop("At least two groups are required for environmental comparisons.",
         call. = FALSE)
  }

  long <- do.call(rbind, lapply(variables, function(v) {
    data.frame(
      SampleID = rownames(env),
      Variable = v,
      Value = as.numeric(env[[v]]),
      Group = grp,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }))
  names(long)[names(long) == "Group"] <- group
  long$Variable <- factor(long$Variable, levels = variables)
  long[[group]] <- factor(long[[group]], levels = levels(grp))
  rownames(long) <- NULL

  overall_list <- list()
  pairwise_list <- list()

  for (v in variables) {
    z <- long[as.character(long$Variable) == v, , drop = FALSE]
    z <- z[is.finite(z$Value) & !is.na(z[[group]]), , drop = FALSE]
    z[[group]] <- droplevels(factor(z[[group]]))
    lev <- levels(z[[group]])
    n_by_group <- table(z[[group]])

    if (length(unique(z$Value)) < 2L) {
      overall_list[[v]] <- data.frame(
        Variable = v, Test = test, N = nrow(z),
        p_value = NA_real_, p_adjusted = NA_real_,
        Note = "Insufficient variation", stringsAsFactors = FALSE
      )
      next
    }
    if (length(lev) < 2L) {
      overall_list[[v]] <- data.frame(
        Variable = v, Test = test, N = nrow(z),
        p_value = NA_real_, p_adjusted = NA_real_,
        Note = "Fewer than two groups with observations", stringsAsFactors = FALSE
      )
      next
    }

    pv <- NA_real_
    note <- ""
    if (test == "kruskal") {
      pv <- tryCatch(stats::kruskal.test(
        stats::reformulate(group, response = "Value"), data = z
      )$p.value, error = function(e) NA_real_)
    } else if (test == "anova") {
      pv <- tryCatch(stats::oneway.test(
        stats::reformulate(group, response = "Value"),
        data = z, var.equal = TRUE
      )$p.value, error = function(e) NA_real_)
    } else {
      note <- "Overall test not requested"
    }

    overall_list[[v]] <- data.frame(
      Variable = v, Test = test, N = nrow(z),
      p_value = as.numeric(pv), p_adjusted = NA_real_,
      Note = note, stringsAsFactors = FALSE
    )

    if (posthoc != "none") {
      pairs <- utils::combn(lev, 2L, simplify = FALSE)
      pp <- lapply(pairs, function(pair) {
        x <- z$Value[z[[group]] == pair[1]]
        y <- z$Value[z[[group]] == pair[2]]
        if (!length(x) || !length(y)) return(NULL)

        dunn_p <- NA_real_
        wilcox_p <- NA_real_

        if (posthoc %in% c("dunn", "both")) {
          all_x <- z$Value
          ranks <- rank(all_x, ties.method = "average")
          n <- length(all_x)
          tie_sizes <- as.numeric(table(all_x))
          tie_correction <- if (n > 1) {
            1 - sum(tie_sizes^3 - tie_sizes) / (n^3 - n)
          } else 0
          if (tie_correction > 0) {
            r1 <- mean(ranks[z[[group]] == pair[1]])
            r2 <- mean(ranks[z[[group]] == pair[2]])
            variance <- n * (n + 1) / 12 * tie_correction *
              (1 / length(x) + 1 / length(y))
            zz <- (r1 - r2) / sqrt(variance)
            dunn_p <- 2 * stats::pnorm(-abs(zz))
          }
        }

        if (posthoc %in% c("wilcox", "both")) {
          wilcox_p <- tryCatch(stats::wilcox.test(
            x, y, exact = FALSE
          )$p.value, error = function(e) NA_real_)
        }

        data.frame(
          Variable = v,
          Group1 = pair[1],
          Group2 = pair[2],
          n1 = length(x),
          n2 = length(y),
          dunn_p = dunn_p,
          wilcox_p = wilcox_p,
          stringsAsFactors = FALSE
        )
      })
      pp <- Filter(Negate(is.null), pp)
      if (length(pp)) {
        pw <- do.call(rbind, pp)
        if (posthoc %in% c("dunn", "both")) {
          pw$dunn_p_adjusted <- stats::p.adjust(pw$dunn_p, method = p_adjust)
        }
        if (posthoc %in% c("wilcox", "both")) {
          pw$wilcox_p_adjusted <- stats::p.adjust(pw$wilcox_p, method = p_adjust)
        }
        pairwise_list[[v]] <- pw
      }
    }
  }

  overall <- do.call(rbind, overall_list)
  rownames(overall) <- NULL
  overall$p_adjusted <- stats::p.adjust(overall$p_value, method = p_adjust)

  pairwise <- if (length(pairwise_list)) {
    out <- do.call(rbind, pairwise_list)
    rownames(out) <- NULL
    out
  } else {
    data.frame()
  }

  structure(
    list(
      data = long,
      statistics = overall,
      pairwise = pairwise,
      parameters = list(
        group = group,
        variables = variables,
        test = test,
        posthoc = posthoc,
        p_adjust = p_adjust,
        group_levels = levels(grp)
      )
    ),
    class = "aq_env"
  )
}

#' Plot environmental factor distributions and group differences
#'
#' @param data An `aq_env` result.
#' @param type `boxplot` or `bar`.
#' @param variables Optional subset of variables stored in `data`.
#' @param facet Whether to facet variables into panels.
#' @param facet_scales Facet scale setting, normally `free_y`.
#' @param bar_error Bar error type: `SE`, `SD`, or `CI` (95%).
#' @param point Show raw sample points.
#' @param point_size Raw sample point size.
#' @param alpha Box/bar transparency.
#' @param palette AqueCo palette.
#' @param colors Optional named group colors.
#' @param show_p Show adjusted overall P values.
#' @param show_sig Show significance stars for overall tests.
#' @param pairwise_label Whether to display adjusted pairwise P values
#'   (or significance stars) above groups.
#' @param pairwise_method Which pairwise result to annotate: `dunn` or `wilcox`.
#' @param pairwise_label_type `p` for adjusted P values or `signif` for stars.
#' @param theme AqueCo theme style.
#' @return An `aq_plot` object with plot and statistics.
#' @export
aq_env_plot <- function(data, type = c("boxplot", "bar"), variables = NULL,
                        facet = TRUE, facet_scales = "free_y",
                        bar_error = c("SE", "SD", "CI"),
                        point = TRUE, point_size = 1.8, alpha = 0.75,
                        palette = "nature", colors = NULL,
                        show_p = TRUE, show_sig = TRUE,
                        pairwise_label = FALSE,
                        pairwise_method = c("dunn", "wilcox"),
                        pairwise_label_type = c("p", "signif"),
                        theme = "publication") {
  if (!inherits(data, "aq_env")) {
    stop("`data` must be an aq_env object.", call. = FALSE)
  }
  type <- match.arg(type)
  bar_error <- match.arg(bar_error)
  pairwise_method <- match.arg(pairwise_method)
  pairwise_label_type <- match.arg(pairwise_label_type)

  d <- data$data
  if (is.null(variables)) variables <- data$parameters$variables
  missing_vars <- setdiff(variables, data$parameters$variables)
  if (length(missing_vars)) {
    stop("Unknown variable(s) in aq_env result: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  d <- d[as.character(d$Variable) %in% variables, , drop = FALSE]
  d$Variable <- factor(as.character(d$Variable), levels = variables)
  group <- data$parameters$group
  d[[group]] <- droplevels(factor(d[[group]]))
  group_levels <- levels(d[[group]])

  cols <- if (is.null(colors)) aq_pal(palette, length(group_levels)) else colors
  if (is.null(names(cols))) names(cols) <- group_levels
  if (!all(group_levels %in% names(cols))) {
    stop("`colors` must be a named vector covering all groups: ",
         paste(group_levels, collapse = ", "), call. = FALSE)
  }
  cols <- cols[group_levels]

  if (type == "boxplot") {
    p <- ggplot2::ggplot(
      d, ggplot2::aes(x = .data[[group]], y = Value,
                      color = .data[[group]])
    ) +
      ggplot2::geom_boxplot(
        ggplot2::aes(fill = .data[[group]]),
        alpha = alpha, outlier.shape = NA, width = 0.65
      )
    if (point) {
      p <- p + ggplot2::geom_jitter(
        width = 0.13, size = point_size, alpha = 0.65
      )
    }
    p <- p + ggplot2::scale_color_manual(values = cols) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(x = group, y = "Environmental value")
  } else {
    summ <- do.call(rbind, lapply(split(d, d$Variable), function(z) {
      spl <- split(z$Value, z[[group]], drop = TRUE)
      do.call(rbind, lapply(names(spl), function(g) {
        x <- spl[[g]]
        x <- x[is.finite(x)]
        n <- length(x)
        m <- if (n) mean(x) else NA_real_
        s <- if (n > 1) stats::sd(x) else NA_real_
        err <- switch(
          bar_error,
          SE = if (n > 1) s / sqrt(n) else NA_real_,
          SD = s,
          CI = if (n > 1) stats::qt(0.975, df = n - 1) * s / sqrt(n)
               else NA_real_
        )
        data.frame(
          Variable = as.character(z$Variable[1]),
          GroupValue = g, Mean = m, Error = err, N = n,
          stringsAsFactors = FALSE
        )
      }))
    }))
    summ$Variable <- factor(summ$Variable, levels = variables)
    summ$GroupValue <- factor(summ$GroupValue, levels = group_levels)
    p <- ggplot2::ggplot(
      summ,
      ggplot2::aes(x = GroupValue, y = Mean, fill = GroupValue)
    ) +
      ggplot2::geom_col(alpha = alpha, width = 0.7, na.rm = TRUE) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data$Mean - .data$Error, ymax = .data$Mean + .data$Error),
        width = 0.18, na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(values = cols) +
      ggplot2::labs(x = group, y = paste0("Mean \u00b1 ", bar_error))
    if (point) {
      p <- p + ggplot2::geom_jitter(
        data = d,
        mapping = ggplot2::aes(
          x = .data[[group]], y = Value,
          color = .data[[group]]
        ),
        inherit.aes = FALSE,
        width = 0.10, size = point_size, alpha = 0.55
      ) +
        ggplot2::scale_color_manual(values = cols, guide = "none")
    }
  }

  if (facet && length(variables) > 1L) {
    p <- p + ggplot2::facet_wrap(
      ~Variable, scales = facet_scales
    )
  } else if (length(variables) == 1L) {
    p <- p + ggplot2::labs(title = variables)
  }

  overall <- data$statistics[
    data$statistics$Variable %in% variables, , drop = FALSE
  ]
  if (show_p && nrow(overall)) {
    overall$label <- if (show_sig) {
      ifelse(is.na(overall$p_adjusted), "NA",
             ifelse(overall$p_adjusted < 0.001, "***",
                    ifelse(overall$p_adjusted < 0.01, "**",
                           ifelse(overall$p_adjusted < 0.05, "*", "ns"))))
    } else {
      paste0("BH-adjusted p = ", signif(overall$p_adjusted, 3))
    }
    ymax <- stats::aggregate(Value ~ Variable, d, max, na.rm = TRUE)
    overall <- merge(overall, ymax, by = "Variable", all.x = TRUE, sort = FALSE)
    overall$y_pos <- overall$Value + 0.08 * pmax(abs(overall$Value), 1)
    overall$LabelX <- factor(utils::tail(group_levels, 1), levels = group_levels)
    p <- p + ggplot2::geom_text(
      data = overall,
      ggplot2::aes(x = LabelX, y = y_pos, label = label),
      inherit.aes = FALSE, hjust = 1, vjust = 0
    )
  }

  if (pairwise_label && nrow(data$pairwise)) {
    pw <- data$pairwise[data$pairwise$Variable %in% variables, , drop = FALSE]
    pcol <- if (pairwise_method == "dunn") "dunn_p_adjusted" else "wilcox_p_adjusted"
    if (!pcol %in% names(pw)) {
      stop("The requested pairwise method was not calculated. Re-run aq_env() ",
           "with posthoc = '", pairwise_method, "' or 'both'.", call. = FALSE)
    }
    pw$P_adjusted <- pw[[pcol]]
    pw$label <- if (pairwise_label_type == "p") {
      paste0("p=", signif(pw$P_adjusted, 3))
    } else {
      ifelse(is.na(pw$P_adjusted), "NA",
             ifelse(pw$P_adjusted < 0.001, "***",
                    ifelse(pw$P_adjusted < 0.01, "**",
                           ifelse(pw$P_adjusted < 0.05, "*", "ns"))))
    }
    # Draw compact pairwise labels at the first group position. This keeps
    # annotations readable in multi-panel plots without crowded bracket stacks.
    if (nrow(pw)) {
      ymax <- stats::aggregate(Value ~ Variable, d, max, na.rm = TRUE)
      pw <- merge(pw, ymax, by = "Variable", all.x = TRUE, sort = FALSE)
      pw <- pw[order(pw$Variable, pw$Group1, pw$Group2), , drop = FALSE]
      pw$pair_index <- stats::ave(seq_len(nrow(pw)), pw$Variable, FUN = seq_along)
      pw$y_pos <- pw$Value +
        (0.08 * pmax(abs(pw$Value), 1)) * pw$pair_index
      pw$Group1 <- factor(pw$Group1, levels = group_levels)
      p <- p + ggplot2::geom_text(
        data = pw,
        ggplot2::aes(x = Group1, y = y_pos, label = label),
        inherit.aes = FALSE, vjust = 0, size = 3, color = "grey20"
      )
    }
  }

  p <- p + theme_aq(theme)
  structure(
    list(
      plot = p,
      statistics = overall,
      pairwise = data$pairwise,
      data = d,
      parameters = list(
        type = type, variables = variables, facet = facet,
        bar_error = if (type == "bar") bar_error else NULL
      )
    ),
    class = "aq_plot"
  )
}

#' Export environmental factor results
#'
#' @param x An `aq_env` object.
#' @param file Output CSV filename prefix.
#' @return File paths, invisibly.
#' @export
aq_env_export <- function(x, file = "environment") {
  if (!inherits(x, "aq_env")) stop("x must be an aq_env object.")
  base <- tools::file_path_sans_ext(file)
  f_data <- paste0(base, "_data.csv")
  f_stats <- paste0(base, "_statistics.csv")
  f_pair <- paste0(base, "_pairwise.csv")
  utils::write.csv(x$data, f_data, row.names = FALSE, na = "")
  utils::write.csv(x$statistics, f_stats, row.names = FALSE, na = "")
  utils::write.csv(x$pairwise, f_pair, row.names = FALSE, na = "")
  invisible(c(data = f_data, statistics = f_stats, pairwise = f_pair))
}

#' @exportS3Method print aq_env
print.aq_env <- function(x, ...) {
  cat(
    "AqueCo environmental factor analysis\n",
    "Grouping variable: ", x$parameters$group, "\n",
    "Variables: ", length(x$parameters$variables), "\n",
    "Overall test: ", x$parameters$test, "\n",
    "Post-hoc: ", x$parameters$posthoc, "\n",
    sep = ""
  )
  print(x$statistics, row.names = FALSE)
  invisible(x)
}
