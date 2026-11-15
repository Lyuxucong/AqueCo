#' Calculate correlations among environmental factors
#'
#' @param data An aq_data object with env data.
#' @param variables Optional numeric environmental variables; NULL selects all.
#' @param method Correlation method: pearson, spearman, or kendall.
#' @param group Optional metadata column for within-group correlations.
#' @param p_adjust P-value adjustment method.
#' @param use Missing-data handling: pairwise.complete.obs or complete.obs.
#' @return An aq_env_cor object containing correlation, P-value, and tidy tables.
#' @export
aq_env_cor <- function(data, variables = NULL,
                       method = c("spearman", "pearson", "kendall"),
                       group = NULL, p_adjust = "BH",
                       use = c("pairwise.complete.obs", "complete.obs")) {
  if (!inherits(data, "aq_data")) stop("`data` must be an aq_data object.", call. = FALSE)
  if (is.null(data$env)) stop("No environmental data found in data$env.", call. = FALSE)
  method <- match.arg(method)
  use <- match.arg(use)
  valid_adjust <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!p_adjust %in% valid_adjust) stop("Invalid `p_adjust` method.", call. = FALSE)

  env <- as.data.frame(data$env, check.names = FALSE)
  if (is.null(rownames(env)) || !all(rownames(data$otu) %in% rownames(env))) {
    stop("data$env row names must contain every sample ID in data$otu.", call. = FALSE)
  }
  env <- env[rownames(data$otu), , drop = FALSE]
  if (is.null(variables)) variables <- names(env)[vapply(env, is.numeric, logical(1))]
  if (!is.character(variables) || length(variables) < 2L) {
    stop("At least two numeric environmental variables must be selected.", call. = FALSE)
  }
  miss <- setdiff(variables, names(env))
  if (length(miss)) stop("Unknown variable(s): ", paste(miss, collapse = ", "), call. = FALSE)
  bad <- variables[!vapply(env[variables], is.numeric, logical(1))]
  if (length(bad)) stop("Variables must be numeric: ", paste(bad, collapse = ", "), call. = FALSE)
  env <- env[, unique(variables), drop = FALSE]

  grp <- NULL
  if (!is.null(group)) {
    if (length(group) != 1L || !group %in% names(data$metadata)) {
      stop("`group` must name one column in data$metadata.", call. = FALSE)
    }
    md <- data$metadata[rownames(env), , drop = FALSE]
    grp <- droplevels(factor(md[[group]]))
    if (anyNA(grp)) stop("Grouping variable contains missing values.", call. = FALSE)
  }

  calc_one <- function(x, label = "All") {
    k <- ncol(x)
    r <- matrix(NA_real_, k, k, dimnames = list(names(x), names(x)))
    p <- r
    n <- matrix(0L, k, k, dimnames = list(names(x), names(x)))
    for (i in seq_len(k)) for (j in seq_len(k)) {
      ok <- is.finite(x[[i]]) & is.finite(x[[j]])
      if (use == "complete.obs") {
        ok <- stats::complete.cases(x)
      }
      xx <- x[[i]][ok]; yy <- x[[j]][ok]
      n[i,j] <- length(xx)
      if (length(xx) >= 3L && length(unique(xx)) > 1L && length(unique(yy)) > 1L) {
        ct <- tryCatch(stats::cor.test(xx, yy, method = method, exact = FALSE),
                       error = function(e) NULL)
        if (!is.null(ct)) {
          r[i,j] <- unname(ct$estimate)
          p[i,j] <- ct$p.value
        }
      } else if (i == j && length(xx) >= 1L) {
        r[i,j] <- 1
        p[i,j] <- 0
      }
    }
    diag(r) <- 1
    diag(p) <- 0
    p_adj <- matrix(stats::p.adjust(p[upper.tri(p)], method = p_adjust),
                    nrow = k, ncol = k)
    p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]
    diag(p_adj) <- 0
    list(r = r, p = p, p_adj = p_adj, n = n, group = label)
  }

  pieces <- if (is.null(grp)) {
    list(All = calc_one(env))
  } else {
    split_idx <- split(seq_len(nrow(env)), grp, drop = TRUE)
    lapply(names(split_idx), function(g) calc_one(env[split_idx[[g]], , drop = FALSE], g)) |>
      stats::setNames(names(split_idx))
  }

  tidy <- do.call(rbind, lapply(names(pieces), function(g) {
    z <- pieces[[g]]
    idx <- which(upper.tri(z$r), arr.ind = TRUE)
    data.frame(
      Group = g,
      Variable1 = rownames(z$r)[idx[,1]],
      Variable2 = colnames(z$r)[idx[,2]],
      Correlation = z$r[idx],
      p_value = z$p[idx],
      p_adjusted = z$p_adj[idx],
      N = z$n[idx],
      stringsAsFactors = FALSE
    )
  }))
  rownames(tidy) <- NULL
  structure(list(
    correlation = lapply(pieces, `[[`, "r"),
    p_value = lapply(pieces, `[[`, "p"),
    p_adjusted = lapply(pieces, `[[`, "p_adj"),
    n = lapply(pieces, `[[`, "n"),
    results = tidy,
    parameters = list(variables = variables, method = method,
                      group = group, p_adjust = p_adjust, use = use)
  ), class = "aq_env_cor")
}

#' Plot environmental-factor correlations as a heatmap
#'
#' @param data An aq_env_cor result.
#' @param show_value Display correlation coefficients in cells.
#' @param show_sig Display significance stars based on adjusted P values.
#' @param label_size Text size in tiles.
#' @param low,mid,high Diverging colors for negative, zero, and positive values.
#' @param theme AqueCo theme.
#' @return An aq_plot object.
#' @export
aq_env_cor_plot <- function(data, show_value = TRUE, show_sig = TRUE,
                            label_size = 3, low = "#3C5488",
                            mid = "#F7F7F7", high = "#E64B35",
                            theme = "publication") {
  if (!inherits(data, "aq_env_cor")) stop("`data` must be an aq_env_cor object.")
  groups <- names(data$correlation)
  tiles <- do.call(rbind, lapply(groups, function(g) {
    r <- data$correlation[[g]]
    pa <- data$p_adjusted[[g]]
    d <- expand.grid(Variable1 = rownames(r), Variable2 = colnames(r),
                     stringsAsFactors = FALSE)
    d$Correlation <- as.vector(r)
    d$p_adjusted <- as.vector(pa)
    d$Group <- g
    d$label <- if (show_value) sprintf("%.2f", d$Correlation) else ""
    if (show_sig) {
      stars <- ifelse(is.na(d$p_adjusted), "",
        ifelse(d$p_adjusted < 0.001, "***",
          ifelse(d$p_adjusted < 0.01, "**",
            ifelse(d$p_adjusted < 0.05, "*", ""))))
      d$label <- paste0(d$label, stars)
    }
    d
  }))
  vars <- data$parameters$variables
  tiles$Variable1 <- factor(tiles$Variable1, levels = vars)
  tiles$Variable2 <- factor(tiles$Variable2, levels = rev(vars))
  p <- ggplot2::ggplot(tiles, ggplot2::aes(x = Variable1, y = Variable2,
                                            fill = Correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient2(
      low = low, mid = mid, high = high, midpoint = 0,
      limits = c(-1, 1), name = "Correlation", na.value = "grey90"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_equal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ))
  if (show_value || show_sig) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label), size = label_size, color = "grey15"
    )
  }
  if (length(groups) > 1L) p <- p + ggplot2::facet_wrap(~Group)
  p <- p + theme_aq(theme)
  structure(list(plot = p, results = data$results,
                 correlation = data$correlation,
                 p_adjusted = data$p_adjusted,
                 parameters = data$parameters),
            class = "aq_plot")
}

#' @exportS3Method print aq_env_cor
print.aq_env_cor <- function(x, ...) {
  cat("AqueCo environmental-factor correlation analysis\n",
      "Method: ", x$parameters$method, "\n",
      "Variables: ", length(x$parameters$variables), "\n",
      "Groups: ", length(x$correlation), "\n", sep = "")
  print(x$results, row.names = FALSE)
  invisible(x)
}
