.aq_distance <- function(otu, method, binary=FALSE) {
  allowed <- c("bray","jaccard","kulczynski","horn","morisita","manhattan",
               "euclidean","canberra","gower","altGower","binomial","chord",
               "cao","mahalanobis")
  if (!method %in% allowed) stop("Unsupported distance. See vegan::vegdist methods.")
  vegan::vegdist(otu, method=method, binary=binary)
}

#' Calculate beta diversity, ordination and PERMANOVA
#' @param data aq_data.
#' @param method Distance method passed to vegan::vegdist.
#' @param ordination NMDS or PCoA.
#' @param group Optional metadata variable for PERMANOVA.
#' @param permutations Number of permutations.
#' @param transform none, hellinger, log, or relative.
#' @param k Number of ordination dimensions.
#' @param trymax NMDS maximum random starts.
#' @param binary Presence/absence distance option.
#' @param formula Optional PERMANOVA formula, e.g. `~ Treatment + Time`.
#' @export
aq_beta <- function(data, method="bray", ordination=c("NMDS","PCoA"),
                    group=NULL, permutations=999, transform=c("none","hellinger","log","relative"),
                    k=2, trymax=100, binary=FALSE, formula=NULL) {
  if (!inherits(data,"aq_data")) stop("data must be aq_data.")
  ordination <- match.arg(ordination); transform <- match.arg(transform)
  x <- data$otu
  if (anyNA(x) || any(x<0) || any(!is.finite(x))) stop("OTU must be finite and non-negative.")
  if (transform=="hellinger") x <- vegan::decostand(x,"hellinger")
  if (transform=="log") x <- log1p(x)
  if (transform=="relative") x <- sweep(x,1,pmax(rowSums(x),1),"/")
  dist <- .aq_distance(x,method,binary)
  if (ordination=="NMDS") {
    ord <- vegan::metaMDS(dist,k=k,trymax=trymax,trace=0,autotransform=FALSE)
    scores <- as.data.frame(vegan::scores(ord,display="sites",choices=seq_len(k)))
    scores$SampleID <- rownames(scores)
    names(scores)[seq_len(k)] <- paste0("Axis",seq_len(k))
    eig <- NULL; variance <- NULL; stress <- ord$stress
  } else {
    cmd <- stats::cmdscale(dist,k=k,eig=TRUE,add=TRUE)
    scores <- as.data.frame(cmd$points)
    scores$SampleID <- rownames(cmd$points)
    names(scores)[seq_len(ncol(cmd$points))] <- paste0("Axis",seq_len(ncol(cmd$points)))
    # Keep the cmdscale result in the common `ord` slot used by the
    # returned aq_beta object; otherwise the PCoA branch leaves `ord`
    # undefined even though the coordinates were calculated successfully.
    ord <- cmd
    eig <- cmd$eig
    pos <- eig[eig > 0]
    variance <- if (length(pos) && sum(pos) > 0) 100 * pos / sum(pos) else numeric(0)
    stress <- NULL
  }
  if (!is.null(data$metadata)) scores <- cbind(scores, data$metadata[match(scores$SampleID,rownames(data$metadata)),,drop=FALSE])
  permanova <- NULL
  if (!is.null(formula) || !is.null(group)) {
    if (is.null(data$metadata)) stop("Metadata is required for PERMANOVA.")
    md <- data$metadata[rownames(x),,drop=FALSE]
    if (is.null(formula)) {
      if (is.null(group) || !group %in% names(md)) stop("Unknown PERMANOVA group.")
      rhs <- paste(group, collapse = " + ")
    } else {
      if (inherits(formula, "formula")) {
        rhs <- paste(deparse(formula[[length(formula)]]), collapse = "")
        rhs <- sub("^~", "", rhs)
      } else {
        stop("formula must be a one-sided formula, e.g. ~ Treatment + Time.")
      }
    }
    perm_formula <- stats::as.formula(paste("dist ~", rhs), env = environment())
    permanova <- vegan::adonis2(perm_formula, data = md,
                                permutations = permutations, by = "margin")
  }
  structure(list(distance=dist, scores=scores, ordination=ord,
                 eigenvalues=eig, variance=variance, stress=stress,
                 permanova=permanova, parameters=list(method=method,
                 ordination=ordination, transform=transform, group=group,
                 permutations=permutations,k=k), source=data),
            class="aq_beta")
}

#' Plot NMDS or PCoA ordination
#' @param data aq_beta result.
#' @param method NMDS or PCoA.
#' @param color Metadata variable mapped to color.
#' @param shape Metadata variable mapped to shape.
#' @param size Point size.
#' @param point_alpha Point transparency.
#' @param palette Palette.
#' @param colors Named custom colors.
#' @param centroid Add group centroids.
#' @param ellipse Add confidence ellipses.
#' @param ellipse_level Ellipse confidence level.
#' @param label Add sample labels.
#' @param marginal Add top and right group-wise boxplots for Axis1 and Axis2 (requires aplot).
#' @param test Statistical test for axis-wise group differences.
#' @param theme Theme style.
#' @export
aq_beta_plot <- function(data, method=NULL, color=NULL, shape=NULL, size=3,
                         point_alpha=.85, palette="nature", colors=NULL,
                         centroid=FALSE, ellipse=FALSE, ellipse_level=.95,
                         label=FALSE, marginal=FALSE, test=c("kruskal","wilcox","none"),
                         theme="publication") {
  if (!inherits(data,"aq_beta")) stop("data must be aq_beta.")
  if (is.null(method)) method <- data$parameters$ordination
  method <- match.arg(method,c("NMDS","PCoA"))
  test <- match.arg(test)
  d <- data$scores

  if (!all(c("Axis1","Axis2") %in% names(d)))
    stop("At least two ordination axes are required.")
  if (!is.null(color) && (!is.character(color) || length(color) != 1L ||
                          !color %in% names(d)))
    stop("Unknown color variable.")
  if (!is.null(shape) && (!is.character(shape) || length(shape) != 1L ||
                          !shape %in% names(d)))
    stop("Unknown shape variable.")
  if (!is.numeric(size) || length(size) != 1L || size <= 0)
    stop("`size` must be one positive number.")
  if (!is.numeric(point_alpha) || length(point_alpha) != 1L ||
      point_alpha < 0 || point_alpha > 1)
    stop("`point_alpha` must be between 0 and 1.")
  if (!is.numeric(ellipse_level) || length(ellipse_level) != 1L ||
      ellipse_level <= 0 || ellipse_level >= 1)
    stop("`ellipse_level` must be between 0 and 1.")

  # Build mappings programmatically. Do not put conditional expressions
  # inside aes(): ggplot2 can interpret their source text as legend labels.
  mapping <- ggplot2::aes(x=Axis1, y=Axis2)
  if (!is.null(color)) mapping$colour <- rlang::quo(.data[[!!color]])
  if (!is.null(shape)) mapping$shape <- rlang::quo(.data[[!!shape]])

  p <- ggplot2::ggplot(d, mapping) +
    ggplot2::geom_point(size=size, alpha=point_alpha)

  if (!is.null(color)) {
    lev <- unique(as.character(d[[color]]))
    cols <- if (is.null(colors)) aq_pal(palette, length(lev)) else colors
    if (is.null(names(cols))) names(cols) <- lev
    p <- p + ggplot2::scale_color_manual(values=cols, name=color)
    if (ellipse) {
      ellipse_mapping <- ggplot2::aes(group=.data[[color]])
      p <- p + ggplot2::stat_ellipse(
        mapping=ellipse_mapping, data=d, inherit.aes=TRUE,
        type="t", level=ellipse_level, linewidth=.6,
        show.legend=FALSE
      )
    }
    if (centroid) {
      centers <- stats::aggregate(
        d[, c("Axis1", "Axis2"), drop=FALSE],
        list(.aq_group=d[[color]]), mean
      )
      names(centers)[names(centers) == ".aq_group"] <- color
      centroid_mapping <- ggplot2::aes(x=Axis1, y=Axis2)
      centroid_mapping$colour <- rlang::quo(.data[[!!color]])
      p <- p + ggplot2::geom_point(
        data=centers, mapping=centroid_mapping,
        inherit.aes=FALSE, shape=4, size=size*2,
        stroke=1.2, show.legend=FALSE
      )
    }
  }

  if (!is.null(shape)) {
    # The legend title is the actual mapped metadata column name.
    p <- p + ggplot2::scale_shape_discrete(name=shape)
  }

  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label=SampleID), nudge_y=.02,
      show.legend=FALSE
    )
  }

  xlab <- if (method=="PCoA" && length(data$variance)>=1)
    sprintf("PCoA1 (%.1f%%)", data$variance[1]) else paste(method, "1")
  ylab <- if (method=="PCoA" && length(data$variance)>=2)
    sprintf("PCoA2 (%.1f%%)", data$variance[2]) else paste(method, "2")
  ttl <- if (method=="NMDS")
    paste0("NMDS (stress = ", format(data$stress, digits=3), ")") else "PCoA"

  # Explicitly set legend labels to metadata names (e.g. "Group").
  # This prevents deparsed conditional expressions from appearing as titles.
  p <- p + ggplot2::labs(
    x=xlab, y=ylab, title=ttl,
    color=if (!is.null(color)) color else NULL,
    shape=if (!is.null(shape)) shape else NULL
  ) + theme_aq(theme)

  axis_tests <- NULL
  if (!is.null(color) && test != "none") {
    axis_tests <- do.call(rbind, lapply(c("Axis1", "Axis2"), function(ax) {
      f <- stats::reformulate(color, response=ax)
      pv <- if (test=="wilcox" && length(unique(d[[color]]))==2)
        stats::wilcox.test(f, data=d)$p.value else
        stats::kruskal.test(f, data=d)$p.value
      data.frame(axis=ax, p_value=pv)
    }))
  }
  if (marginal) {
    if (!requireNamespace("aplot", quietly=TRUE)) {
      stop(
        "Marginal boxplots require the 'aplot' package. ",
        "Install it with install.packages('aplot').",
        call. = FALSE
      )
    }
    if (is.null(color)) {
      stop(
        "Marginal boxplots require `color` to specify the grouping variable, ",
        "for example color = \"Group\".",
        call. = FALSE
      )
    }

    # Keep group order identical to the main ordination plot.
    group_values <- d[[color]]
    if (!is.factor(group_values)) {
      group_values <- factor(
        as.character(group_values),
        levels = unique(as.character(group_values))
      )
    }
    marginal_data <- d
    marginal_data$.aq_marginal_group <- group_values

    group_levels <- levels(group_values)
    marginal_cols <- if (!is.null(colors)) colors else
      aq_pal(palette, length(group_levels))
    if (is.null(names(marginal_cols))) {
      names(marginal_cols) <- group_levels
    }
    marginal_cols <- marginal_cols[group_levels]

    # Top panel: horizontal boxplots of Axis1. coord_flip() places the
    # numerical Axis1 scale on the horizontal direction so it aligns with
    # the main ordination x-axis, following the microeco/aplot layout.
    p_top <- ggplot2::ggplot(
      marginal_data,
      ggplot2::aes(
        x = .aq_marginal_group, y = Axis1,
        fill = .aq_marginal_group, colour = .aq_marginal_group
      )
    ) +
      ggplot2::geom_boxplot(
        width = 0.62, alpha = 0.55, linewidth = 0.35,
        outlier.shape = NA
      ) +
      ggplot2::scale_fill_manual(values = marginal_cols, guide = "none") +
      ggplot2::scale_colour_manual(values = marginal_cols, guide = "none") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(2, 5.5, 0, 5.5)
      )

    # Right panel: vertical boxplots of Axis2, aligned to the main y-axis.
    p_right <- ggplot2::ggplot(
      marginal_data,
      ggplot2::aes(
        x = .aq_marginal_group, y = Axis2,
        fill = .aq_marginal_group, colour = .aq_marginal_group
      )
    ) +
      ggplot2::geom_boxplot(
        width = 0.62, alpha = 0.55, linewidth = 0.35,
        outlier.shape = NA
      ) +
      ggplot2::scale_fill_manual(values = marginal_cols, guide = "none") +
      ggplot2::scale_colour_manual(values = marginal_cols, guide = "none") +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5.5, 2, 5.5, 0)
      )

    # Remove the title from the base ggplot BEFORE aplot composes it.
    # Adding ggplot layers after insert_top()/insert_right() can dispatch
    # to a composite-object method and may discard the composed plot.
    p <- p + ggplot2::labs(title = NULL) +
      ggplot2::theme(plot.title = ggplot2::element_blank())

    # Compose the marginal panels last and retain each returned object.
    p <- aplot::insert_top(p, p_top, height = 0.20)
    p <- aplot::insert_right(p, p_right, width = 0.20)

    if (is.null(p)) {
      stop(
        "The marginal plot composition returned NULL. Check the installed ",
        "version of 'aplot' and its insert_top()/insert_right() methods.",
        call. = FALSE
      )
    }
  }
  structure(
    list(plot=p, axis_tests=axis_tests, permanova=data$permanova,
         scores=d, stress=data$stress, variance=data$variance),
    class="aq_plot"
  )
}

#' @exportS3Method print aq_beta
print.aq_beta <- function(x,...) {
  cat("AqueCo beta diversity:",x$parameters$ordination,
      "\nDistance:",x$parameters$method,
      if (!is.null(x$stress)) paste0("\nStress: ",signif(x$stress,4)) else "",
      "\n")
  if (!is.null(x$permanova)) print(x$permanova)
  invisible(x)
}
