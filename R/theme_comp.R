#' AqueCo ggplot theme
#' @param style Theme style. `publication` and `bw` use theme_bw().
#' @param base_size Base font size.
#' @param legend_position Legend position.
#' @export
theme_aq <- function(style = c("publication", "bw", "classic", "minimal", "nature"),
                     base_size = 12, legend_position = "right") {
  style <- match.arg(style)
  th <- switch(style,
    publication = ggplot2::theme_bw(base_size = base_size),
    bw = ggplot2::theme_bw(base_size = base_size),
    classic = ggplot2::theme_classic(base_size = base_size),
    minimal = ggplot2::theme_minimal(base_size = base_size),
    nature = ggplot2::theme_bw(base_size = base_size)
  )
  th + ggplot2::theme(
    legend.position = legend_position,
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold"),
    # Keep AqueCo's standard ggplot themes compatible with ggtern.
    # ggtern registers this special theme element; a non-rel value inherited
    # from the active global theme can otherwise break ordinary ggplot plots.
    tern.axis.ticks.length.major = ggplot2::rel(0.5),
    tern.axis.ticks.length.minor = ggplot2::rel(0.25)
  )
}

#' Set the global AqueCo plotting theme
#' @param style Theme style passed to theme_aq().
#' @param base_size Base font size.
#' @param legend_position Legend position.
#' @return The previous ggplot2 theme, invisibly.
#' @export
aq_theme_set <- function(style = "publication", base_size = 12,
                         legend_position = "right") {
  ggplot2::theme_set(
    theme_aq(style = style, base_size = base_size,
             legend_position = legend_position)
  )
}

#' Plot taxonomic composition
#' @param data aq_data or aq_taxa.
#' @param type Plot type: bar, heatmap, bubble, or venn.
#' @param level Taxonomic rank when data is aq_data.
#' @param summary sample or mean.
#' @param group Metadata grouping variable; required for summary = "mean".
#' @param facet Metadata faceting variable(s).
#' @param abundance count or relative.
#' @param top Optional top N taxa ranked by mean relative abundance.
#' @param other Combine taxa outside top N into Other (bar only).
#' @param palette Palette name.
#' @param colors Optional named discrete colors (bar/bubble).
#' @param position stack or fill (bar only).
#' @param theme Theme style.
#' @param facet_scales Facet scale setting.
#' @param x_angle Angle of x-axis labels in degrees.
#' @param x_size X-axis label text size.
#' @param x_labels Whether to display x-axis labels.
#' @param heatmap_scale Heatmap transformation: none or row_z.
#' @param heatmap_colors Three colors for low, middle, and high heatmap values.
#' @param cluster_samples Cluster samples (ignored for summary = "mean").
#' @param cluster_taxa Cluster taxa.
#' @param bubble_color Bubble color mapping: taxon, group, or abundance.
#' @param bubble_size_range Minimum and maximum bubble sizes.
#' @param alpha Point/tile transparency for bubble plots.
#' @param show_values Display values in heatmap tiles.
#' @param prevalence Minimum within-group prevalence for Venn sets.
#' @param min_abundance Minimum abundance for Venn set membership.
#' @export
aq_comp <- function(data, type = c("bar", "heatmap", "bubble", "venn"),
                    level = "Genus",
                    summary = c("sample", "mean"), group = NULL, facet = NULL,
                    abundance = "relative", top = NULL, other = FALSE,
                    palette = "nature", colors = NULL,
                    position = c("stack", "fill"), theme = "publication",
                    facet_scales = c("free_x", "fixed", "free", "free_y"),
                    x_angle = 90, x_size = 8, x_labels = TRUE,
                    heatmap_scale = c("none", "row_z"),
                    heatmap_colors = c("#3C5488", "#F7F7F7", "#E64B35"),
                    cluster_samples = FALSE, cluster_taxa = FALSE,
                    bubble_color = c("taxon", "group", "abundance"),
                    bubble_size_range = c(1.5, 9),
                    alpha = 0.85, show_values = FALSE,
                    prevalence = 0, min_abundance = 0) {

  type <- match.arg(type)
  summary <- match.arg(summary)
  position <- match.arg(position)
  facet_scales <- match.arg(facet_scales)
  heatmap_scale <- match.arg(heatmap_scale)
  bubble_color <- match.arg(bubble_color)

  # Venn is a set-comparison plot and returns the richer aq_venn result object.
  if (type == "venn") {
    if (is.null(group) || length(group) != 1L)
      stop("type = 'venn' requires group to name one metadata variable.")
    return(aq_venn(
      data = data, group = group, level = level,
      prevalence = prevalence, min_abundance = min_abundance,
      palette = palette
    ))
  }

  if (!is.numeric(x_angle) || length(x_angle) != 1 ||
      !is.finite(x_angle) || x_angle < 0 || x_angle > 180)
    stop("x_angle must be between 0 and 180.")
  if (!is.numeric(x_size) || length(x_size) != 1 ||
      !is.finite(x_size) || x_size <= 0)
    stop("x_size must be a positive number.")
  if (!is.logical(x_labels) || length(x_labels) != 1 || is.na(x_labels))
    stop("x_labels must be TRUE or FALSE.")
  if (!is.logical(cluster_samples) || length(cluster_samples) != 1 ||
      is.na(cluster_samples) ||
      !is.logical(cluster_taxa) || length(cluster_taxa) != 1 ||
      is.na(cluster_taxa))
    stop("cluster_samples and cluster_taxa must be TRUE or FALSE.")
  if (!is.numeric(alpha) || length(alpha) != 1 ||
      !is.finite(alpha) || alpha < 0 || alpha > 1)
    stop("alpha must be between 0 and 1.")
  if (!is.numeric(bubble_size_range) || length(bubble_size_range) != 2 ||
      any(!is.finite(bubble_size_range)) ||
      bubble_size_range[1] < 0 || bubble_size_range[2] <= bubble_size_range[1])
    stop("bubble_size_range must contain two increasing non-negative values.")
  if (length(heatmap_colors) != 3)
    stop("heatmap_colors must contain exactly three colors: low, middle, high.")

  tx <- if (inherits(data, "aq_data")) {
    aq_taxa(data, level = level, abundance = "count")
  } else {
    data
  }
  if (!inherits(tx, "aq_taxa"))
    stop("data must be an aq_data or aq_taxa object.")

  abundance <- match.arg(abundance, c("count", "relative"))
  mat <- if (abundance == "count") tx$count else tx$relative
  if (is.null(rownames(mat)) || is.null(colnames(mat)) ||
      !nrow(mat) || !ncol(mat))
    stop("The abundance matrix is empty or lacks dimnames.")

  # Rank taxa by mean relative abundance, regardless of displayed abundance.
  if (!is.null(top)) {
    if (length(top) != 1 || !is.numeric(top) || is.na(top) ||
        top < 1 || top != as.integer(top))
      stop("top must be a positive integer.")
    ranking <- colMeans(tx$relative, na.rm = TRUE)
    keep <- names(sort(ranking, decreasing = TRUE))[
      seq_len(min(as.integer(top), length(ranking)))
    ]
    mat <- mat[, intersect(keep, colnames(mat)), drop = FALSE]
    if (type == "bar" && isTRUE(other)) {
      drop_taxa <- setdiff(colnames(if (abundance == "count") tx$count else tx$relative),
                           colnames(mat))
      source_mat <- if (abundance == "count") tx$count else tx$relative
      if (length(drop_taxa)) {
        mat <- cbind(
          mat,
          Other = rowSums(source_mat[, drop_taxa, drop = FALSE], na.rm = TRUE)
        )
      }
    }
  }

  if (!is.null(tx$metadata)) {
    missing_meta <- setdiff(unique(c(group, facet)), names(tx$metadata))
    missing_meta <- missing_meta[!is.na(missing_meta) & nzchar(missing_meta)]
    if (length(missing_meta))
      stop("Unknown metadata variable(s): ", paste(missing_meta, collapse = ", "))
  } else if (!is.null(group) || !is.null(facet)) {
    stop("Metadata is required when group or facet is specified.")
  }

  # Optional clustering is performed on the sample x taxon matrix.
  if (summary == "sample" && isTRUE(cluster_samples) && nrow(mat) > 1) {
    sample_order <- stats::hclust(stats::dist(mat))$order
    mat <- mat[sample_order, , drop = FALSE]
  }
  if (isTRUE(cluster_taxa) && ncol(mat) > 1) {
    tax_order <- stats::hclust(stats::dist(t(mat)))$order
    mat <- mat[, tax_order, drop = FALSE]
  }

  d <- expand.grid(
    SampleID = rownames(mat),
    Taxon = colnames(mat),
    stringsAsFactors = FALSE
  )
  d$Abundance <- as.vector(mat)

  meta <- tx$metadata
  if (!is.null(meta)) {
    meta <- meta[match(d$SampleID, rownames(meta)), , drop = FALSE]
    d <- cbind(d, meta)
  }

  if (summary == "mean") {
    if (is.null(group) || length(group) != 1 || !group %in% names(d))
      stop("Specify one valid group for summary = 'mean'.")
    by <- unique(c(group, facet, "Taxon"))
    d <- stats::aggregate(
      Abundance ~ .,
      data = d[, c(by, "Abundance"), drop = FALSE],
      FUN = mean
    )
    x <- group
  } else {
    x <- "SampleID"
  }

  # Keep x-axis order stable and preserve the order of samples/groups.
  if (summary == "sample") {
    d$SampleID <- factor(d$SampleID, levels = rownames(mat))
  } else {
    d[[x]] <- factor(d[[x]], levels = unique(as.character(d[[x]])))
  }
  d$Taxon <- factor(d$Taxon, levels = colnames(mat))

  if (type == "bar") {
    cols <- if (is.null(colors)) {
      stats::setNames(aq_pal(palette, n = nlevels(d$Taxon)), levels(d$Taxon))
    } else colors

    p <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = .data[[x]], y = Abundance, fill = Taxon)
    ) +
      ggplot2::geom_col(position = position) +
      ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
      ggplot2::labs(
        x = x, y = abundance,
        fill = if (!is.null(tx$parameters$level)) tx$parameters$level else level
      ) +
      theme_aq(theme)

  } else if (type == "heatmap") {
    if (heatmap_scale == "row_z") {
      # Standardize each taxon across the displayed x-axis categories.
      # This also works after summary = "mean", where d has fewer rows
      # than the original sample-by-taxon matrix.
      d$PlotValue <- stats::ave(d$Abundance, d$Taxon, FUN = function(z) {
        s <- stats::sd(z, na.rm = TRUE)
        if (!is.finite(s) || s == 0) rep(0, length(z)) else
          (z - mean(z, na.rm = TRUE)) / s
      })
      fill_label <- "Row Z-score"
    } else {
      d$PlotValue <- d$Abundance
      fill_label <- if (abundance == "relative") "Relative abundance" else "Count"
    }

    p <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = .data[[x]], y = Taxon, fill = PlotValue)
    ) +
      ggplot2::geom_tile(
        color = "white",
        linewidth = 0.2,
        alpha = alpha
      ) +
      ggplot2::scale_fill_gradient2(
        low = heatmap_colors[1],
        mid = heatmap_colors[2],
        high = heatmap_colors[3],
        midpoint = if (heatmap_scale == "row_z") 0 else
          stats::median(d$PlotValue, na.rm = TRUE),
        name = fill_label
      ) +
      ggplot2::labs(x = x, y = if (!is.null(tx$parameters$level)) tx$parameters$level else level) +
      theme_aq(theme)

    if (isTRUE(show_values)) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = formatC(PlotValue, digits = 2, format = "f")),
        size = 2.5,
        color = "#263238"
      )
    }

  } else {
    if (bubble_color == "group" &&
        (is.null(group) || !group %in% names(d)))
      stop("bubble_color = 'group' requires a valid group variable.")
    if (bubble_color == "taxon") {
      color_mapping <- ggplot2::aes(color = Taxon)
      bubble_scale <- ggplot2::scale_color_manual(
        values = if (is.null(colors)) {
          stats::setNames(aq_pal(palette, n = nlevels(d$Taxon)), levels(d$Taxon))
        } else colors
      )
    } else if (bubble_color == "group") {
      group_values <- unique(as.character(d[[group]]))
      color_mapping <- ggplot2::aes(color = .data[[group]])
      bubble_scale <- ggplot2::scale_color_manual(
        values = if (is.null(colors)) {
          stats::setNames(aq_pal(palette, n = length(group_values)), group_values)
        } else colors
      )
    } else {
      color_mapping <- ggplot2::aes(color = Abundance)
      bubble_scale <- ggplot2::scale_color_gradientn(
        colors = aq_pal(palette, n = 100, type = "continuous")
      )
    }

    p <- ggplot2::ggplot(
      d,
      ggplot2::aes(x = .data[[x]], y = Taxon, size = Abundance)
    ) +
      ggplot2::geom_point(
        mapping = color_mapping,
        alpha = alpha
      ) +
      ggplot2::scale_size_continuous(
        range = bubble_size_range,
        name = if (abundance == "relative") "Relative abundance" else "Count"
      ) +
      bubble_scale +
      ggplot2::labs(x = x, y = if (!is.null(tx$parameters$level)) tx$parameters$level else level) +
      theme_aq(theme)
  }

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(
      stats::reformulate(facet),
      scales = facet_scales,
      drop = TRUE
    )
  }

  p <- p + ggplot2::theme(
    axis.text.x = if (x_labels) {
      ggplot2::element_text(
        angle = x_angle,
        hjust = if (x_angle == 0) 0.5 else 1,
        vjust = 0.5,
        size = x_size
      )
    } else {
      ggplot2::element_blank()
    }
  )
  p
}
