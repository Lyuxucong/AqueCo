.aq_sets <- function(data, group, level = "Genus", prevalence = 0, min_abundance = 0) {
  tx <- if (inherits(data, "aq_data")) aq_taxa(data, level, "count") else data
  if (!inherits(tx, "aq_taxa")) stop("data must be aq_data or aq_taxa.")
  meta <- tx$metadata
  if (is.null(meta) || !group %in% names(meta)) stop("A valid metadata group is required.")
  g <- droplevels(factor(meta[[group]]))
  if (anyNA(g) || nlevels(g) < 2) stop("At least two non-missing groups are required.")
  sets <- lapply(levels(g), function(z) {
    m <- tx$count[g == z, , drop = FALSE] > min_abundance
    names(which(colMeans(m) >= prevalence & colSums(m) > 0))
  })
  names(sets) <- levels(g)
  sets
}

#' Create a Venn diagram
#' @param data aq_data or aq_taxa.
#' @param group Metadata grouping variable.
#' @param level Taxonomic level.
#' @param prevalence Within-group prevalence threshold.
#' @param min_abundance Minimum count to call presence.
#' @param palette Palette.
#' @export
aq_venn <- function(data, group, level = "Genus", prevalence = 0,
                    min_abundance = 0, palette = "nature") {
  if (!requireNamespace("ggVennDiagram", quietly = TRUE))
    stop("Install ggVennDiagram to use aq_venn().")
  if (prevalence < 0 || prevalence > 1) stop("prevalence must be in [0,1].")
  sets <- .aq_sets(data, group, level, prevalence, min_abundance)
  if (length(sets) > 5) stop("Use aq_upset() for more than five groups.")
  p <- ggVennDiagram::ggVennDiagram(sets) +
    ggplot2::scale_fill_gradientn(colours = aq_pal(palette, 64, "continuous")) +
    theme_aq()
  ints <- list()
  nms <- names(sets)
  for (k in seq_along(nms)) for (cmb in utils::combn(nms, k, simplify = FALSE)) {
    z <- Reduce(intersect, sets[cmb])
    ex <- setdiff(nms, cmb)
    if (length(ex)) z <- setdiff(z, unlist(sets[ex]))
    ints[[paste(cmb, collapse = " & ")]] <- z
  }
  uniq <- lapply(nms, function(g) setdiff(sets[[g]], unlist(sets[setdiff(nms, g)])))
  names(uniq) <- nms
  structure(list(plot = p, sets = sets, unique = uniq,
                 shared = Reduce(intersect, sets), intersections = ints,
                 intersection_table = data.frame(intersection = names(ints), n_taxa = lengths(ints)),
                 parameters = list(group = group, level = level, prevalence = prevalence)),
            class = "aq_venn")
}

#' Create an UpSet plot
#' @param data aq_data or aq_taxa.
#' @param group Metadata grouping variable.
#' @param level Taxonomic rank.
#' @param prevalence Within-group prevalence.
#' @param min_abundance Minimum abundance.
#' @param top Number of intersections.
#' @export
aq_upset <- function(data, group, level = "Genus", prevalence = 0,
                     min_abundance = 0, top = 20) {
  if (!requireNamespace("ComplexUpset", quietly = TRUE))
    stop("Install ComplexUpset to use aq_upset().")
  sets <- .aq_sets(data, group, level, prevalence, min_abundance)
  taxa <- unique(unlist(sets))
  d <- data.frame(Taxon = taxa, check.names = FALSE)
  for (nm in names(sets)) d[[nm]] <- taxa %in% sets[[nm]]
  p <- ComplexUpset::upset(d, intersect = names(sets), n_intersections = top)
  structure(list(plot = p, sets = sets), class = "aq_upset")
}

#' Export taxa intersections
#' @param x aq_venn or aq_upset result.
#' @param file CSV path.
#' @export
aq_sets_export <- function(x, file = "aq_intersections.csv") {
  if (is.null(x$intersections)) stop("This result does not contain intersection memberships.")
  d <- do.call(rbind, lapply(names(x$intersections), function(nm) {
    z <- x$intersections[[nm]]
    if (!length(z)) return(NULL)
    data.frame(Intersection = nm, Taxon = z)
  }))
  if (is.null(d)) d <- data.frame(Intersection = character(), Taxon = character())
  utils::write.csv(d, file, row.names = FALSE)
  invisible(file)
}
