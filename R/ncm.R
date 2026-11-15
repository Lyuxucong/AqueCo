
#' Fit the Sloan neutral community model
#'
#' Fits a separate Sloan neutral community model to the whole community or
#' independently within each metadata group.
#' @param data An aq_data object.
#' @param group Optional metadata column defining independent communities.
#' @param groups Optional subset of levels to analyze.
#' @param level Optional taxonomic rank to aggregate with aq_taxa.
#' @param min_prevalence Minimum within-community prevalence for retaining taxa.
#' @param min_abundance Minimum count defining taxon presence.
#' @param N Optional community size used in the Sloan model. Defaults to mean
#'   library size within each analysis set.
#' @return An aq_ncm object containing summary, taxon-level results and fits.
#' @export
aq_ncm <- function(data, group = NULL, groups = NULL, level = NULL,
                   min_prevalence = 0, min_abundance = 0, N = NULL) {
  if (!inherits(data, "aq_data")) stop("`data` must be an aq_data object.", call. = FALSE)
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1L ||
      !is.finite(min_prevalence) || min_prevalence < 0 || min_prevalence > 1)
    stop("`min_prevalence` must be between 0 and 1.", call. = FALSE)
  if (!is.numeric(min_abundance) || length(min_abundance) != 1L ||
      !is.finite(min_abundance) || min_abundance < 0)
    stop("`min_abundance` must be a non-negative number.", call. = FALSE)
  comm <- if (is.null(level)) as.matrix(data$otu) else
    as.matrix(aq_taxa(data, level = level, abundance = "count")$abundance)
  if (any(!is.finite(comm)) || any(comm < 0))
    stop("Abundance matrix must contain finite non-negative values.", call. = FALSE)
  if (is.null(group)) {
    if (!is.null(groups)) stop("`groups` can only be used when `group` is supplied.")
    md_group <- rep("All", nrow(comm))
    names(md_group) <- rownames(comm)
    group_levels <- "All"
  } else {
    if (is.null(data$metadata)) stop("Metadata are required when `group` is specified.")
    if (length(group) != 1L || !group %in% names(data$metadata))
      stop("`group` must name a column in data$metadata.")
    md_group <- as.character(data$metadata[rownames(comm), group])
    if (anyNA(md_group) || any(!nzchar(md_group))) stop("Grouping variable contains missing/empty values.")
    group_levels <- unique(md_group)
    if (!is.null(groups)) {
      unknown <- setdiff(as.character(groups), group_levels)
      if (length(unknown)) stop("Unknown group(s): ", paste(unknown, collapse = ", "))
      group_levels <- as.character(groups)
    }
  }
  if (!is.null(N) && (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0))
    stop("`N` must be one positive number or NULL.")
  fits <- list(); taxa_out <- list(); summaries <- list()
  for (g in group_levels) {
    idx <- which(md_group == g)
    x <- comm[idx, , drop = FALSE]
    x <- x[rowSums(x) > 0, , drop = FALSE]
    if (nrow(x) < 3L) {
      warning("Group ", g, " has fewer than 3 non-empty samples; skipped.")
      next
    }
    lib <- rowSums(x)
    n <- nrow(x)
    prev <- colMeans(x > min_abundance)
    p <- colSums(x) / sum(x)
    keep <- is.finite(prev) & prev >= min_prevalence & p > 0 & p < 1
    prev <- prev[keep]; p <- p[keep]
    if (length(p) < 2L) {
      warning("Group ", g, " has fewer than 2 eligible taxa; skipped.")
      next
    }
    Ng <- if (is.null(N)) mean(lib) else N
    pred <- function(log_nm) {
      nm <- exp(log_nm)
      1 - stats::pbeta(1 / Ng, nm * p, nm * (1 - p))
    }
    objective <- function(log_nm) {
      yhat <- pred(log_nm)
      sum((prev - yhat)^2)
    }
    opt <- stats::optimize(objective, interval = log(c(1e-6, 1e9)))
    nm <- exp(opt$minimum)
    fitted <- pred(opt$minimum)
    r2 <- if (stats::var(prev) > 0) 1 - sum((prev - fitted)^2) /
      sum((prev - mean(prev))^2) else NA_real_
    se <- sqrt(pmax(fitted * (1 - fitted) / n, 0))
    lo <- pmax(0, fitted - stats::qnorm(0.975) * se)
    hi <- pmin(1, fitted + stats::qnorm(0.975) * se)
    cls <- ifelse(prev > hi, "Above", ifelse(prev < lo, "Below", "Neutral"))
    td <- data.frame(Group = g, Taxon = names(p), MeanRelativeAbundance = unname(p),
      OccurrenceFrequency = unname(prev), PredictedFrequency = unname(fitted),
      CI_lower = unname(lo), CI_upper = unname(hi), Classification = cls,
      stringsAsFactors = FALSE)
    fits[[g]] <- list(Nm = nm, N = Ng, m = nm / Ng, R2 = r2,
                      n_samples = n, n_taxa = length(p), data = td)
    taxa_out[[g]] <- td
    summaries[[g]] <- data.frame(Group = g, R2 = r2, Nm = nm, m = nm / Ng,
                                  N = Ng, Samples = n, Taxa = length(p))
  }
  if (!length(fits)) stop("No groups could be fitted; check sample numbers and filtering.")
  structure(list(summary = do.call(rbind, summaries),
    taxa = do.call(rbind, taxa_out), fits = fits,
    parameters = list(group = group, groups = group_levels, level = level,
      min_prevalence = min_prevalence, min_abundance = min_abundance,
      N = N, method = "Sloan neutral community model",
      interval = "Approximate 95% binomial interval")), class = "aq_ncm")
}

#' Plot Sloan neutral community model results
#' @param data An aq_ncm result.
#' @param facet Whether to facet groups in one plot.
#' @param show_ci Show approximate 95% binomial confidence limits.
#' @param show_class Color taxa by Above, Neutral and Below classes.
#' @param point_size Point size.
#' @param point_alpha Point transparency.
#' @param palette Named colors for Above, Neutral and Below.
#' @param theme AqueCo theme style.
#' @return A ggplot object, or a named list of ggplot objects when facet is FALSE.
#' @export
aq_ncm_plot <- function(data, facet = TRUE, show_ci = TRUE, show_class = TRUE,
                        point_size = 2, point_alpha = 0.65,
                        palette = c(Above = "#D95F5F", Neutral = "#8A8A8A",
                                    Below = "#4C91C6"), theme = "publication") {
  if (!inherits(data, "aq_ncm")) stop("`data` must be an aq_ncm object.")
  make_plot <- function(d, g) {
    s <- data$summary[data$summary$Group == g, , drop = FALSE]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = MeanRelativeAbundance,
      y = OccurrenceFrequency)) +
      ggplot2::geom_line(ggplot2::aes(y = PredictedFrequency),
        color = "#245B78", linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Classification),
        size = point_size, alpha = point_alpha) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = c(0, .02))) +
      ggplot2::labs(title = g, x = "Mean relative abundance",
        y = "Occurrence frequency", color = "Model deviation") +
      ggplot2::annotate("text", x = Inf, y = Inf,
        label = sprintf("R\u00b2 = %.3f\nNm = %.2f", s$R2, s$Nm),
        hjust = 1.08, vjust = 1.2, size = 3.8) +
      theme_aq(theme) +
      ggplot2::theme(legend.position = "bottom")
    if (isTRUE(show_ci)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
        fill = "#9ECAE1", alpha = 0.35, inherit.aes = TRUE,
        na.rm = TRUE
      )
    }
    if (!isTRUE(show_class)) p <- p + ggplot2::guides(color = "none")
    if (!is.null(palette)) p <- p + ggplot2::scale_color_manual(values = palette, drop = FALSE)
    if (!isTRUE(show_ci)) p <- p + ggplot2::geom_blank()
    p
  }
  if (isTRUE(facet)) {
    d <- data$taxa
    p <- ggplot2::ggplot(d, ggplot2::aes(x = MeanRelativeAbundance,
      y = OccurrenceFrequency)) +
      ggplot2::geom_line(ggplot2::aes(y = PredictedFrequency), color = "#245B78",
        linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(color = Classification),
        size = point_size, alpha = point_alpha) +
      ggplot2::facet_wrap(~Group, scales = "fixed") +
      ggplot2::scale_x_log10() + ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = "Mean relative abundance", y = "Occurrence frequency",
        color = "Model deviation") + theme_aq(theme) +
      ggplot2::theme(legend.position = "bottom")
    if (isTRUE(show_ci)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = CI_lower, ymax = CI_upper),
        fill = "#9ECAE1", alpha = 0.3
      )
    }
    ann <- data$summary
    ann$label <- sprintf("R\u00b2 = %.3f\nNm = %.2f", ann$R2, ann$Nm)
    ann$x <- Inf; ann$y <- Inf
    p <- p + ggplot2::geom_text(data = ann,
      ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE,
      hjust = 1.08, vjust = 1.2, size = 3.5)
    if (!is.null(palette)) p <- p + ggplot2::scale_color_manual(values = palette, drop = FALSE)
    return(p)
  }
  stats::setNames(lapply(names(data$fits), function(g)
    make_plot(data$fits[[g]]$data, g)), names(data$fits))
}

#' Export NCM results
#' @param data An aq_ncm result.
#' @param file Output file prefix.
#' @return Invisibly returns the output file paths.
#' @export
aq_ncm_export <- function(data, file = "NCM_results") {
  if (!inherits(data, "aq_ncm")) stop("`data` must be an aq_ncm object.")
  summary_file <- paste0(file, "_summary.csv")
  taxa_file <- paste0(file, "_taxa.csv")
  utils::write.csv(data$summary, summary_file, row.names = FALSE)
  utils::write.csv(data$taxa, taxa_file, row.names = FALSE)
  invisible(c(summary = summary_file, taxa = taxa_file))
}
