#' Correlate environmental factors with microbial taxa
#'
#' Taxa are aggregated at the requested taxonomy level and correlated with
#' numeric environmental variables across matched samples.
#'
#' @param data An aq_data object containing env and taxonomy.
#' @param env_variables Optional environmental columns; NULL selects all numeric.
#' @param level Taxonomic rank used to aggregate OTUs/ASVs.
#' @param abundance Taxon abundance scale: relative or count.
#' @param method Pearson, Spearman, or Kendall correlation.
#' @param group Optional metadata column; correlations are calculated separately
#'   within each group.
#' @param prevalence Minimum fraction of samples in the analysis subset in which
#'   a taxon must be present (0 to 1).
#' @param min_abundance Minimum abundance used to define presence.
#' @param top Optional maximum number of taxa, ranked by mean abundance.
#' @param p_adjust P-value adjustment method across all tested taxon-environment
#'   pairs within each analysis group.
#' @param use Missing-data handling: pairwise.complete.obs or complete.obs.
#' @return An aq_env_microbe_cor object with pairwise results and abundance data.
#' @export
aq_env_microbe_cor <- function(data, env_variables = NULL, level = "Genus",
                               abundance = c("relative", "count"),
                               method = c("spearman", "pearson", "kendall"),
                               group = NULL, prevalence = 0.10,
                               min_abundance = 0, top = NULL,
                               p_adjust = "BH",
                               use = c("pairwise.complete.obs", "complete.obs")) {
  if (!inherits(data, "aq_data")) stop("`data` must be an aq_data object.", call. = FALSE)
  if (is.null(data$env)) stop("Environmental data are required in data$env.", call. = FALSE)
  if (is.null(data$taxonomy)) stop("Taxonomy is required to aggregate microbial taxa.", call. = FALSE)
  abundance <- match.arg(abundance)
  method <- match.arg(method)
  use <- match.arg(use)
  valid_adjust <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!is.character(p_adjust) || length(p_adjust) != 1L || !p_adjust %in% valid_adjust)
    stop("Invalid `p_adjust` method.", call. = FALSE)
  if (!is.numeric(prevalence) || length(prevalence) != 1L ||
      !is.finite(prevalence) || prevalence < 0 || prevalence > 1)
    stop("`prevalence` must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(min_abundance) || length(min_abundance) != 1L ||
      !is.finite(min_abundance) || min_abundance < 0)
    stop("`min_abundance` must be a non-negative number.", call. = FALSE)

  env <- as.data.frame(data$env, check.names = FALSE)
  if (is.null(rownames(env)) || !all(rownames(data$otu) %in% rownames(env)))
    stop("data$env row names must contain every sample ID in data$otu.", call. = FALSE)
  env <- env[rownames(data$otu), , drop = FALSE]
  if (is.null(env_variables)) {
    env_variables <- names(env)[vapply(env, is.numeric, logical(1))]
  }
  if (!is.character(env_variables) || !length(env_variables))
    stop("At least one numeric environmental variable must be selected.", call. = FALSE)
  miss_env <- setdiff(env_variables, names(env))
  if (length(miss_env)) stop("Unknown environmental variable(s): ",
                             paste(miss_env, collapse = ", "), call. = FALSE)
  nonnum <- env_variables[!vapply(env[env_variables], is.numeric, logical(1))]
  if (length(nonnum)) stop("Environmental variables must be numeric: ",
                           paste(nonnum, collapse = ", "), call. = FALSE)
  env <- env[, unique(env_variables), drop = FALSE]

  tax_obj <- aq_taxa(data, level = level, abundance = abundance, unclassified = "keep")
  taxa_abund <- tax_obj$abundance
  taxa_abund <- taxa_abund[rownames(env), , drop = FALSE]

  grp <- NULL
  if (!is.null(group)) {
    if (is.null(data$metadata) || length(group) != 1L ||
        !group %in% names(data$metadata))
      stop("`group` must name a column in data$metadata.", call. = FALSE)
    md <- data$metadata[rownames(env), , drop = FALSE]
    grp <- droplevels(factor(md[[group]]))
    if (anyNA(grp)) stop("Grouping variable contains missing values.", call. = FALSE)
  } else {
    grp <- factor(rep("All", nrow(env)), levels = "All")
  }

  result_list <- list()
  abundance_list <- list()
  for (g in levels(grp)) {
    idx <- which(grp == g)
    eg <- env[idx, , drop = FALSE]
    tg <- taxa_abund[idx, , drop = FALSE]
    prev <- colMeans(tg > min_abundance, na.rm = TRUE)
    keep <- is.finite(prev) & prev >= prevalence
    tg <- tg[, keep, drop = FALSE]
    if (!ncol(tg)) next
    if (!is.null(top)) {
      if (!is.numeric(top) || length(top) != 1L || top < 1)
        stop("`top` must be NULL or a positive integer.", call. = FALSE)
      ord <- order(colMeans(tg, na.rm = TRUE), decreasing = TRUE)
      tg <- tg[, utils::head(ord, as.integer(top)), drop = FALSE]
    }
    abundance_list[[g]] <- data.frame(
      SampleID = rownames(tg),
      Group = g,
      tg,
      check.names = FALSE
    )
    pairs <- expand.grid(
      Taxon = colnames(tg), Environment = colnames(eg),
      stringsAsFactors = FALSE
    )
    one <- lapply(seq_len(nrow(pairs)), function(k) {
      tx <- pairs$Taxon[k]; ev <- pairs$Environment[k]
      x <- as.numeric(tg[, tx]); y <- as.numeric(eg[[ev]])
      ok <- if (use == "complete.obs") stats::complete.cases(tg, eg) else
        is.finite(x) & is.finite(y)
      x <- x[ok]; y <- y[ok]
      n <- length(x)
      rho <- pv <- NA_real_
      if (n >= 3L && length(unique(x)) > 1L && length(unique(y)) > 1L) {
        ct <- tryCatch(stats::cor.test(x, y, method = method, exact = FALSE),
                       error = function(e) NULL)
        if (!is.null(ct)) {
          rho <- unname(ct$estimate)
          pv <- ct$p.value
        }
      }
      data.frame(
        Group = g, Taxon = tx, Environment = ev,
        Correlation = rho, p_value = pv, N = n,
        stringsAsFactors = FALSE
      )
    })
    rr <- do.call(rbind, one)
    rr$p_adjusted <- stats::p.adjust(rr$p_value, method = p_adjust)
    result_list[[g]] <- rr
  }
  if (!length(result_list))
    stop("No taxa passed the prevalence filter. Reduce `prevalence` or check abundance data.",
         call. = FALSE)
  results <- do.call(rbind, result_list)
  rownames(results) <- NULL
  # Taxa can pass prevalence filtering in one group but not another.
  # Align group-specific abundance tables to the union of taxon columns
  # before row-binding; NA means that taxon was not retained for that group.
  all_taxa <- unique(unlist(lapply(abundance_list, function(x) {
    setdiff(names(x), c("SampleID", "Group"))
  }), use.names = FALSE))
  abundance_list <- lapply(abundance_list, function(x) {
    missing_taxa <- setdiff(all_taxa, names(x))
    for (tx in missing_taxa) x[[tx]] <- NA_real_
    x[, c("SampleID", "Group", all_taxa), drop = FALSE]
  })
  abund_out <- do.call(rbind, abundance_list)
  rownames(abund_out) <- NULL
  structure(list(
    results = results,
    abundance = abund_out,
    parameters = list(
      env_variables = env_variables, level = level, abundance = abundance,
      method = method, group = group, prevalence = prevalence,
      min_abundance = min_abundance, top = top,
      p_adjust = p_adjust, use = use
    )
  ), class = "aq_env_microbe_cor")
}

#' Plot environmental factor and microbial taxon correlations
#'
#' @param data An aq_env_microbe_cor result.
#' @param show_value Display correlation coefficients.
#' @param show_sig Display adjusted-P significance stars.
#' @param only_significant Display only adjusted-P significant pairs.
#' @param p_cutoff Adjusted-P cutoff when only_significant is TRUE.
#' @param label_size Text size for cell labels.
#' @param low,mid,high Colors for negative, zero, and positive correlations.
#' @param theme AqueCo theme.
#' @return An aq_plot object.
#' @export
aq_env_microbe_cor_plot <- function(data, show_value = TRUE, show_sig = TRUE,
                                    only_significant = FALSE, p_cutoff = 0.05,
                                    label_size = 3, low = "#3C5488",
                                    mid = "#F7F7F7", high = "#E64B35",
                                    theme = "publication") {
  if (!inherits(data, "aq_env_microbe_cor"))
    stop("`data` must be an aq_env_microbe_cor object.", call. = FALSE)
  d <- data$results
  if (only_significant) d <- d[!is.na(d$p_adjusted) & d$p_adjusted < p_cutoff, , drop = FALSE]
  if (!nrow(d)) stop("No correlations to plot after filtering.", call. = FALSE)
  d$Taxon <- factor(d$Taxon, levels = rev(unique(data$results$Taxon)))
  d$Environment <- factor(d$Environment, levels = unique(data$results$Environment))
  d$label <- if (show_value) sprintf("%.2f", d$Correlation) else ""
  if (show_sig) {
    stars <- ifelse(is.na(d$p_adjusted), "",
      ifelse(d$p_adjusted < 0.001, "***",
        ifelse(d$p_adjusted < 0.01, "**",
          ifelse(d$p_adjusted < 0.05, "*", ""))))
    d$label <- paste0(d$label, stars)
  }
  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = Environment, y = Taxon,
                    fill = Correlation, size = abs(Correlation))
  ) +
    ggplot2::geom_point(shape = 21, color = "grey30", alpha = 0.9) +
    ggplot2::scale_fill_gradient2(
      low = low, mid = mid, high = high, midpoint = 0,
      limits = c(-1, 1), name = "Correlation", na.value = "grey90"
    ) +
    ggplot2::scale_size_continuous(range = c(2, 8), limits = c(0, 1),
                                   name = "|Correlation|") +
    ggplot2::labs(x = "Environmental factors", y = data$parameters$level) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    ))
  if (show_value || show_sig) {
    p <- p + ggplot2::geom_text(
      data = d,
      mapping = ggplot2::aes(x = Environment, y = Taxon, label = label),
      inherit.aes = FALSE, size = label_size, color = "grey15"
    )
  }
  if (length(unique(d$Group)) > 1L || !is.null(data$parameters$group))
    p <- p + ggplot2::facet_wrap(~Group)
  p <- p + theme_aq(theme)
  structure(list(plot = p, results = data$results,
                 data = d, parameters = data$parameters),
            class = "aq_plot")
}

#' @exportS3Method print aq_env_microbe_cor
print.aq_env_microbe_cor <- function(x, ...) {
  cat("AqueCo environmental factor\u2013microbial taxon correlation\n",
      "Taxonomic level: ", x$parameters$level, "\n",
      "Method: ", x$parameters$method, "\n",
      "Tested associations: ", nrow(x$results), "\n", sep = "")
  print(utils::head(x$results), row.names = FALSE)
  invisible(x)
}
