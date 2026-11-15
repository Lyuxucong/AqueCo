#' Perform group-wise Mantel tests between microbial communities and environment
#'
#' For each group, computes a microbial community dissimilarity matrix and
#' tests it against the Euclidean distance matrix of each selected environmental
#' variable. Also calculates the Pearson correlation matrix among environmental
#' variables for the companion lower-triangle plot.
#'
#' @param data An aq_data object containing otu, metadata, and env.
#' @param group Metadata column defining treatment groups.
#' @param env_variables Environmental variables to test; NULL selects all numeric.
#' @param level Taxonomic rank for aggregating microbial abundances; NULL uses OTUs.
#' @param abundance Relative abundance or counts when level is not NULL.
#' @param community_method Community dissimilarity, e.g. bray, jaccard, euclidean.
#' @param env_method Distance for each environmental variable: euclidean or manhattan.
#' @param permutations Number of permutations.
#' @param p_adjust P-value adjustment method within each group.
#' @param prevalence Minimum taxon prevalence within each group.
#' @param min_abundance Minimum abundance defining presence.
#' @param cor_method Correlation method for the environmental matrix.
#' @return An aq_mantel object with Mantel results and environmental correlations.
#' @export
aq_mantel <- function(data, group, env_variables = NULL,
                      level = NULL, abundance = c("relative", "count"),
                      community_method = "bray", env_method = "euclidean",
                      permutations = 999, p_adjust = "BH",
                      prevalence = 0.10, min_abundance = 0,
                      cor_method = "pearson") {
  if (!inherits(data, "aq_data")) stop("`data` must be an aq_data object.", call. = FALSE)
  if (is.null(data$env)) stop("Environmental data are required in data$env.", call. = FALSE)
  if (is.null(data$metadata)) stop("Metadata are required for group-wise Mantel tests.", call. = FALSE)
  if (missing(group) || length(group) != 1L || !group %in% names(data$metadata))
    stop("`group` must name one column in data$metadata.", call. = FALSE)
  if (!is.numeric(permutations) || length(permutations) != 1L ||
      permutations < 1 || permutations != as.integer(permutations))
    stop("`permutations` must be a positive integer.", call. = FALSE)
  if (!is.numeric(prevalence) || length(prevalence) != 1L ||
      prevalence < 0 || prevalence > 1)
    stop("`prevalence` must be between 0 and 1.", call. = FALSE)
  abundance <- match.arg(abundance)
  valid_adjust <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!p_adjust %in% valid_adjust) stop("Invalid p_adjust method.", call. = FALSE)

  env <- as.data.frame(data$env, check.names = FALSE)
  samples <- rownames(data$otu)
  if (is.null(rownames(env)) || !all(samples %in% rownames(env)))
    stop("data$env row names must contain all sample IDs.", call. = FALSE)
  env <- env[samples, , drop = FALSE]
  if (is.null(env_variables))
    env_variables <- names(env)[vapply(env, is.numeric, logical(1))]
  if (!is.character(env_variables) || length(env_variables) < 1L)
    stop("Select at least one numeric environmental variable.", call. = FALSE)
  unknown <- setdiff(env_variables, names(env))
  if (length(unknown)) stop("Unknown environmental variables: ",
                            paste(unknown, collapse = ", "), call. = FALSE)
  bad <- env_variables[!vapply(env[env_variables], is.numeric, logical(1))]
  if (length(bad)) stop("Environmental variables must be numeric: ",
                        paste(bad, collapse = ", "), call. = FALSE)
  env <- env[, unique(env_variables), drop = FALSE]
  if (ncol(env) < 2L)
    warning("Only one environmental variable selected; Pearson matrix will be 1 x 1.",
            call. = FALSE)

  md <- data$metadata[samples, , drop = FALSE]
  grp <- droplevels(factor(md[[group]]))
  if (anyNA(grp)) stop("Grouping variable contains missing values.", call. = FALSE)
  if (is.null(level)) {
    comm <- as.matrix(data$otu)
    level_used <- "OTU"
  } else {
    tx <- aq_taxa(data, level = level, abundance = abundance)
    comm <- as.matrix(tx$abundance)
    level_used <- tx$parameters$level
  }
  if (any(!is.finite(comm)) || any(comm < 0))
    stop("Microbial abundance matrix contains invalid values.", call. = FALSE)

  # Environmental Pearson correlation matrix is calculated over all matched samples.
  cor_mat <- stats::cor(env, method = cor_method, use = "pairwise.complete.obs")
  cor_p <- matrix(NA_real_, ncol(env), ncol(env),
                  dimnames = list(names(env), names(env)))
  for (i in seq_len(ncol(env))) for (j in seq_len(ncol(env))) {
    ok <- stats::complete.cases(env[[i]], env[[j]])
    if (i == j) {
      cor_p[i, j] <- 0
    } else if (sum(ok) >= 3L && length(unique(env[[i]][ok])) > 1L &&
               length(unique(env[[j]][ok])) > 1L) {
      cor_p[i, j] <- tryCatch(
        stats::cor.test(env[[i]][ok], env[[j]][ok], method = cor_method)$p.value,
        error = function(e) NA_real_
      )
    }
  }
  cor_p_adj <- matrix(NA_real_, ncol(env), ncol(env),
                      dimnames = dimnames(cor_p))
  upper <- upper.tri(cor_p)
  cor_p_adj[upper] <- stats::p.adjust(cor_p[upper], method = p_adjust)
  cor_p_adj[lower.tri(cor_p_adj)] <- t(cor_p_adj)[lower.tri(cor_p_adj)]
  diag(cor_p_adj) <- 0

  res_list <- list()
  group_levels <- levels(grp)
  for (g in group_levels) {
    idx <- which(grp == g)
    if (length(idx) < 3L) {
      warning("Group ", g, " has fewer than 3 samples; Mantel tests skipped.")
      next
    }
    eg <- env[idx, , drop = FALSE]
    cg <- comm[idx, , drop = FALSE]
    prev <- colMeans(cg > min_abundance, na.rm = TRUE)
    cg <- cg[, is.finite(prev) & prev >= prevalence, drop = FALSE]
    if (!ncol(cg)) {
      warning("No taxa passed prevalence filtering in group ", g, ".")
      next
    }
    # Remove taxa with zero variance to avoid degenerate distance matrices.
    cg <- cg[, apply(cg, 2, function(x) length(unique(x)) > 1L), drop = FALSE]
    if (!ncol(cg)) {
      warning("No variable taxa remain in group ", g, ".")
      next
    }
    comm_dist <- tryCatch(
      vegan::vegdist(cg, method = community_method),
      error = function(e) stop("Community distance failed in group ", g, ": ", conditionMessage(e))
    )
    for (v in names(eg)) {
      x <- eg[[v]]
      ok <- is.finite(x)
      if (sum(ok) < 3L || length(unique(x[ok])) < 2L) {
        res_list[[length(res_list) + 1L]] <- data.frame(
          Group = g, Environment = v, Mantel_r = NA_real_,
          Mantel_p = NA_real_, N = sum(ok), Note = "Insufficient variation or samples"
        )
        next
      }
      # Keep the same sample set for both distance matrices.
      keep <- which(ok)
      comm_sub <- cg[keep, , drop = FALSE]
      if (nrow(comm_sub) < 3L) next
      cd <- vegan::vegdist(comm_sub, method = community_method)
      ed <- stats::dist(matrix(x[keep], ncol = 1), method = env_method)
      mt <- tryCatch(
        vegan::mantel(cd, ed, method = "pearson", permutations = as.integer(permutations)),
        error = function(e) NULL
      )
      res_list[[length(res_list) + 1L]] <- data.frame(
        Group = g, Environment = v,
        Mantel_r = if (is.null(mt)) NA_real_ else unname(mt$statistic),
        Mantel_p = if (is.null(mt)) NA_real_ else mt$signif,
        N = length(keep),
        Note = if (is.null(mt)) "Mantel test failed" else "",
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(res_list))
    stop("No group had enough samples and taxa for Mantel testing.", call. = FALSE)
  results <- do.call(rbind, res_list)
  rownames(results) <- NULL
  results$Mantel_p_adjusted <- stats::ave(
    results$Mantel_p, results$Group,
    FUN = function(x) stats::p.adjust(x, method = p_adjust)
  )
  structure(list(
    results = results,
    env_data = env,
    env_correlation = cor_mat,
    env_p_value = cor_p,
    env_p_adjusted = cor_p_adj,
    community = comm,
    parameters = list(
      group = group, group_levels = group_levels,
      env_variables = names(env), level = level_used,
      abundance = if (is.null(level)) "input OTU abundance" else abundance,
      community_method = community_method, env_method = env_method,
      permutations = as.integer(permutations), p_adjust = p_adjust,
      prevalence = prevalence, min_abundance = min_abundance,
      cor_method = cor_method
    )
  ), class = "aq_mantel")
}

#' Plot a combined environmental correlation and group-wise Mantel network
#'
#' @param data An aq_mantel result.
#' @param cor_show_sig Show significance stars in the environmental correlation matrix.
#' @param cor_p_cutoff Adjusted P threshold for correlation stars.
#' @param mantel_p_breaks Two cutoffs separating Mantel P-value categories.
#' @param mantel_p_labels Labels for the three Mantel P categories.
#' @param mantel_r_breaks Cutoffs for Mantel line-width categories.
#' @param show_mantel_label Show Mantel r labels beside links.
#' @param show_group_labels Show treatment-group names at the external link nodes.
#' @param group_label_size Text size for treatment-group node labels.
#' @param low,mid,high Colors for negative, zero, and positive Pearson correlations.
#' @param theme AqueCo theme.
#' @param x_label_angle Angle of the bottom heatmap axis labels, default 90 degrees.
#' @param x_label_size Font size of the bottom heatmap axis labels.
#' @return An aq_plot object.
#' @export
aq_mantel_plot <- function(data, cor_show_sig = TRUE, cor_p_cutoff = 0.05,
                           mantel_p_breaks = c(0.01, 0.05),
                           mantel_p_labels = c("< 0.01", "0.01 \u2013 0.05", ">= 0.05"),
                           mantel_r_breaks = c(0.2, 0.3),
                           show_mantel_label = FALSE,
                           show_group_labels = TRUE,
                           group_label_size = 3.8,
                           low = "#3C80B8", mid = "#F7F7F7",
                           high = "#F04438", theme = "publication",
                           x_label_angle = 90, x_label_size = 10) {
  if (!inherits(data, "aq_mantel")) stop("`data` must be an aq_mantel object.")
  if (!is.numeric(x_label_angle) || length(x_label_angle) != 1L ||
      !is.finite(x_label_angle) || x_label_angle < 0 || x_label_angle > 360)
    stop("`x_label_angle` must be a finite number between 0 and 360.")
  if (!is.numeric(x_label_size) || length(x_label_size) != 1L ||
      !is.finite(x_label_size) || x_label_size <= 0)
    stop("`x_label_size` must be a positive finite number.")
  if (!requireNamespace("linkET", quietly = TRUE)) {
    stop(
      "Package `linkET` is required for Mantel network plotting. Install it with: ",
      "remotes::install_github('Hy4m/linkET')",
      call. = FALSE
    )
  }
  env <- as.data.frame(data$env_data, check.names = FALSE)
  if (is.null(env) || ncol(env) < 2L)
    stop("At least two environmental variables are required for the correlation matrix.")
  if (length(mantel_p_breaks) != 2L || is.unsorted(mantel_p_breaks))
    stop("`mantel_p_breaks` must contain two increasing values.")
  if (length(mantel_r_breaks) != 2L || is.unsorted(mantel_r_breaks))
    stop("`mantel_r_breaks` must contain two increasing values.")

  # linkET geom_couple connects external treatment groups (from)
  # to environmental variables in the matrix (to), without manual node coordinates.
  link_dat <- data$results
  link_dat$from <- as.character(link_dat$Group)
  link_dat$to <- as.character(link_dat$Environment)
  link_dat$p_category <- cut(
    link_dat$Mantel_p, breaks = c(-Inf, mantel_p_breaks, Inf),
    labels = mantel_p_labels, right = FALSE
  )
  link_dat$r_category <- cut(
    abs(link_dat$Mantel_r), breaks = c(-Inf, mantel_r_breaks, Inf),
    labels = c("< 0.2", "0.2 \u2013 0.3", "> 0.3"), right = FALSE
  )
  link_dat <- link_dat[
    !is.na(link_dat$Mantel_r) & !is.na(link_dat$Mantel_p) &
      !is.na(link_dat$from) & !is.na(link_dat$to), , drop = FALSE
  ]
  if (!nrow(link_dat)) stop("No valid Mantel results are available for plotting.")

  p <- linkET::qcorrplot(
    linkET::correlate(env, method = data$parameters$cor_method),
    type = "lower", diag = FALSE
  ) +
    linkET::geom_square(color = "grey65", size = 0.35) +
    ggplot2::scale_fill_gradient2(
      low = low, mid = mid, high = high, midpoint = 0,
      limits = c(-1, 1), name = "Pearson's r", na.value = "grey90"
    )

  if (isTRUE(cor_show_sig)) {
    p <- p + linkET::geom_mark(
      sig_level = c(0.05, 0.01, 0.001),
      sig_thres = cor_p_cutoff, only_mark = TRUE,
      size = 4, color = "black"
    )
  }

  p <- p +
    linkET::geom_couple(
      data = link_dat,
      mapping = ggplot2::aes(
        from = from, to = to,
        colour = p_category, size = r_category
      ),
      curvature = linkET::nice_curvature(),
      label.size = if (isTRUE(show_group_labels)) group_label_size else 0,
      label.colour = "black",
      label.fontface = 1,
      nudge_x = 0.2,
      drop = TRUE
    ) +
    ggplot2::scale_size_manual(
      values = c("< 0.2" = 0.6, "0.2 \u2013 0.3" = 1.2, "> 0.3" = 2.0),
      name = "Mantel's r", drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(
        c("#62C6A5", "#FF805A", "#91A5D6"), mantel_p_labels
      ),
      name = "Mantel's p", drop = FALSE
    ) +
    ggplot2::guides(
      size = ggplot2::guide_legend(
        title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2
      ),
      colour = ggplot2::guide_legend(
        title = "Mantel's p", override.aes = list(size = 1.2), order = 3
      ),
      fill = ggplot2::guide_colorbar(title = "Pearson's r", order = 1)
    ) +
    theme_aq(theme) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = x_label_angle,
        hjust = if (x_label_angle == 0) 0.5 else 1,
        vjust = if (x_label_angle == 0) 0.5 else 0.5,
        size = x_label_size
      )
    )

  structure(
    list(plot = p, results = data$results,
         env_correlation = data$env_correlation,
         env_p_adjusted = data$env_p_adjusted,
         parameters = data$parameters),
    class = "aq_plot"
  )
}

#' @exportS3Method print aq_mantel
print.aq_mantel <- function(x, ...) {
  cat("AqueCo group-wise Mantel analysis\n",
      "Grouping variable: ", x$parameters$group, "\n",
      "Community distance: ", x$parameters$community_method, "\n",
      "Groups: ", length(unique(x$results$Group)), "\n",
      "Environmental variables: ", length(x$parameters$env_variables), "\n",
      sep = "")
  print(x$results, row.names = FALSE)
  invisible(x)
}
