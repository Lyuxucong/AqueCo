#' Decompose beta diversity within metadata groups
#'
#' Computes pairwise Sorensen beta-diversity partitioning independently
#' within each group. Each group contributes its own set of sample pairs,
#' allowing group-specific ternary diagrams to be compared.
#'
#' For each pair, the three ternary components are:
#' turnover = beta.sim; nestedness = beta.sne; similarity = 1 - beta.sor.
#' Their percentages sum to 100.
#'
#' @param data An `aq_data` object (samples in rows, taxa in columns).
#' @param group Metadata column defining groups. Required.
#' @param method Partitioning family: `sorensen` or `jaccard`.
#' @param prevalence Optional within-group prevalence threshold between 0 and 1 (inclusive).
#'   A taxon is retained in a group when its prevalence is >= this value.
#'   The default 0 retains any taxon present in at least one sample.
#' @return An `aq_beta_decomp` object containing pairwise components,
#'   group means, and settings.
#' @export
aq_beta_decomp <- function(data, group, method = c("sorensen", "jaccard"),
                           prevalence = 0) {
  method <- match.arg(method)

  if (!inherits(data, "aq_data")) {
    stop("`data` must be an aq_data object.", call. = FALSE)
  }
  if (missing(group) || !is.character(group) || length(group) != 1L ||
      is.na(group) || !nzchar(group)) {
    stop("Please supply one metadata grouping column through `group`.",
         call. = FALSE)
  }
  if (is.null(data$metadata) || !group %in% names(data$metadata)) {
    stop("Grouping variable `", group, "` was not found in data$metadata.",
         call. = FALSE)
  }
  if (!requireNamespace("betapart", quietly = TRUE)) {
    stop("Install 'betapart' to calculate beta diversity decomposition.",
         call. = FALSE)
  }
  if (!is.numeric(prevalence) || length(prevalence) != 1L ||
      !is.finite(prevalence) || prevalence < 0 || prevalence > 1) {
    stop("`prevalence` must be a single number between 0 and 1.",
         call. = FALSE)
  }

  otu <- as.matrix(data$otu)
  storage.mode(otu) <- "numeric"
  if (anyNA(otu) || any(!is.finite(otu)) || any(otu < 0)) {
    stop("OTU abundance must contain finite, non-negative values.",
         call. = FALSE)
  }

  md <- data$metadata[rownames(otu), , drop = FALSE]
  grp <- md[[group]]
  if (anyNA(grp)) {
    stop("Grouping variable `", group, "` contains missing values.",
         call. = FALSE)
  }
  grp <- droplevels(factor(grp))
  if (nlevels(grp) < 2L) {
    stop("At least two groups are required for group-wise comparison.",
         call. = FALSE)
  }

  # Presence/absence conversion follows the supplied reference workflow.
  pa <- (otu > 0) * 1L
  if (any(rowSums(pa) == 0)) {
    empty <- rownames(pa)[rowSums(pa) == 0]
    stop(
      "The following samples contain no observed taxa and cannot be used: ",
      paste(empty, collapse = ", "),
      call. = FALSE
    )
  }

  group_levels <- levels(grp)
  too_small <- group_levels[vapply(group_levels, function(g) sum(grp == g) < 2L,
                                   logical(1))]
  if (length(too_small)) {
    stop(
      "Each group needs at least two samples for within-group pairwise ",
      "decomposition. Too few samples in: ",
      paste(too_small, collapse = ", "),
      call. = FALSE
    )
  }

  pairwise_list <- vector("list", length(group_levels))
  names(pairwise_list) <- group_levels

  for (g in group_levels) {
    sample_ids <- rownames(pa)[grp == g]
    xg <- pa[sample_ids, , drop = FALSE]

    # Optional prevalence filtering is applied independently in each group.
    if (prevalence > 0) {
      keep_taxa <- colMeans(xg > 0) >= prevalence
      xg <- xg[, keep_taxa, drop = FALSE]
    }
    # betapart requires at least one taxon column; retain a clear error.
    if (ncol(xg) == 0L) {
      stop(
        "No taxa remain in group `", g,
        "` after applying prevalence = ", prevalence, ".",
        call. = FALSE
      )
    }

    bp <- betapart::beta.pair(xg, index.family = method)

    if (method == "sorensen") {
      turnover_mat <- as.matrix(bp$beta.sim)
      nested_mat <- as.matrix(bp$beta.sne)
      total_mat <- as.matrix(bp$beta.sor)
    } else {
      turnover_mat <- as.matrix(bp$beta.jtu)
      nested_mat <- as.matrix(bp$beta.jne)
      total_mat <- as.matrix(bp$beta.jac)
    }

    ij <- which(upper.tri(total_mat), arr.ind = TRUE)
    total <- total_mat[ij]
    turnover <- turnover_mat[ij]
    nestedness <- nested_mat[ij]

    # Similarity is the complement of total dissimilarity. Components
    # therefore form a closed composition (sum = 100%) for ternary plotting.
    similarity <- 1 - total
    out <- data.frame(
      Group = g,
      Sample1 = rownames(total_mat)[ij[, 1]],
      Sample2 = colnames(total_mat)[ij[, 2]],
      Turnover = as.numeric(turnover),
      Nestedness = as.numeric(nestedness),
      TotalDissimilarity = as.numeric(total),
      Similarity = as.numeric(similarity),
      Turnover_pct = as.numeric(turnover) * 100,
      Nestedness_pct = as.numeric(nestedness) * 100,
      Similarity_pct = as.numeric(similarity) * 100,
      stringsAsFactors = FALSE
    )
    pairwise_list[[g]] <- out
  }

  pairwise <- do.call(rbind, pairwise_list)
  rownames(pairwise) <- NULL
  pairwise$Group <- factor(pairwise$Group, levels = group_levels)

  means <- stats::aggregate(
    pairwise[, c("Similarity_pct", "Turnover_pct", "Nestedness_pct")],
    by = list(Group = pairwise$Group),
    FUN = mean,
    na.rm = TRUE
  )

  structure(
    list(
      data = pairwise,
      means = means,
      groups = group_levels,
      parameters = list(
        group = group,
        method = method,
        prevalence = prevalence,
        comparison = "within_group",
        input = "presence_absence"
      ),
      source = data
    ),
    class = "aq_beta_decomp"
  )
}

#' Plot group-wise beta diversity decomposition as ternary diagrams
#'
#' Produces a faceted ternary diagram with one panel per metadata group.
#' Also returns a separate ternary plot for each group in `group_plots`.
#'
#' @param data An `aq_beta_decomp` result.
#' @param point_size Pairwise-comparison point size.
#' @param point_alpha Point transparency between 0 and 1.
#' @param show_mean Whether to show the group mean as a larger diamond.
#' @param mean_size Mean-point size.
#' @param palette AqueCo discrete palette.
#' @param colors Optional named colors for groups.
#' @param theme Theme style passed to `theme_aq()`.
#' @return An `aq_plot` object. The combined faceted plot is in `$plot`;
#'   individual group plots are in `$group_plots`.
#' @export
aq_beta_decomp_plot <- function(data, point_size = 1.8,
                                point_alpha = 0.30, show_mean = TRUE,
                                mean_size = 4, palette = "nature",
                                colors = NULL, theme = "publication") {
  if (!inherits(data, "aq_beta_decomp")) {
    stop("`data` must be an aq_beta_decomp object.", call. = FALSE)
  }
  if (!requireNamespace("ggtern", quietly = TRUE)) {
    stop("Install 'ggtern' to draw ternary diagrams.", call. = FALSE)
  }
  if (!is.numeric(point_size) || length(point_size) != 1L ||
      !is.finite(point_size) || point_size <= 0) {
    stop("`point_size` must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(point_alpha) || length(point_alpha) != 1L ||
      !is.finite(point_alpha) || point_alpha < 0 || point_alpha > 1) {
    stop("`point_alpha` must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(mean_size) || length(mean_size) != 1L ||
      !is.finite(mean_size) || mean_size <= 0) {
    stop("`mean_size` must be a positive number.", call. = FALSE)
  }

  d <- data$data
  means <- data$means
  group_levels <- data$groups

  cols <- if (is.null(colors)) aq_pal(palette, length(group_levels)) else colors
  if (is.null(names(cols))) names(cols) <- group_levels
  if (!all(group_levels %in% names(cols))) {
    stop("When `colors` is supplied, it must be a named vector covering: ",
         paste(group_levels, collapse = ", "), call. = FALSE)
  }
  cols <- cols[group_levels]

  base_plot <- function(dat, mean_dat = NULL, title = NULL) {
    p <- ggtern::ggtern(
      dat,
      ggplot2::aes(
        x = Similarity_pct,
        y = Turnover_pct,
        z = Nestedness_pct
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(color = Group),
        size = point_size,
        alpha = point_alpha
      ) +
      ggplot2::scale_color_manual(values = cols, drop = FALSE) +
      ggtern::Tlab("Similarity (%)") +
      ggtern::Llab("Turnover (%)") +
      ggtern::Rlab("Nestedness (%)") +
      ggplot2::labs(
        title = title,
        color = data$parameters$group
      ) +
      ggtern::theme_bw() +
      ggplot2::theme(
        # ggtern requires this special theme element to be a rel object.
        tern.axis.ticks.length.major = ggplot2::rel(0.5),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        strip.text = ggplot2::element_text(face = "bold"),
        legend.position = "right",
        panel.grid.major = ggplot2::element_line(
          color = "grey82", linetype = "dashed", linewidth = 0.35
        ),
        panel.background = ggplot2::element_rect(fill = "white")
      )

    if (show_mean && !is.null(mean_dat) && nrow(mean_dat)) {
      # ggtern::ggtern() supplies ternary coordinate handling, while the
      # point geom itself is exported by ggplot2 (not ggtern). Inherit the
      # global z mapping and override only x/y for the mean-point data.
      p <- p + ggplot2::geom_point(
  data = mean_dat,
  ggplot2::aes(
    x = Similarity_pct,
    y = Turnover_pct,
    color = Group,
    fill = Group
  ),
  inherit.aes = TRUE,
  shape = 23,
  size = mean_size,
  stroke = 1.1,
  show.legend = FALSE
) +
ggplot2::scale_fill_manual(
  values = cols,
  guide = "none"
)
    }
    p
  }

  # One combined comparison figure: each group gets its own ternary panel.
  p <- base_plot(d, means, title = NULL) +
    ggplot2::facet_wrap(~Group, nrow = 1, drop = FALSE)

  # Also provide separately addressable plots for each group.
  group_plots <- lapply(group_levels, function(g) {
    dg <- d[as.character(d$Group) == g, , drop = FALSE]
    mg <- means[as.character(means$Group) == g, , drop = FALSE]
    base_plot(dg, mg, title = g)
  })
  names(group_plots) <- group_levels

  structure(
    list(
      plot = p,
      group_plots = group_plots,
      data = d,
      means = means,
      parameters = data$parameters
    ),
    class = "aq_plot"
  )
}

#' @exportS3Method print aq_beta_decomp
print.aq_beta_decomp <- function(x, ...) {
  cat(
    "AqueCo group-wise beta diversity decomposition\n",
    "Grouping variable: ", x$parameters$group, "\n",
    "Index family: ", x$parameters$method, "\n",
    "Groups: ", paste(x$groups, collapse = ", "), "\n",
    "Pairwise comparisons: ", nrow(x$data), "\n",
    sep = ""
  )
  print(utils::head(x$data))
  invisible(x)
}
