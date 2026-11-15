.aq_rank_alias <- list(
  Kingdom = c("kingdom", "domain", "k", "d"),
  Phylum = c("phylum", "p"), Class = c("class", "c"),
  Order = c("order", "o"), Family = c("family", "f"),
  Genus = c("genus", "g"), Species = c("species", "s")
)

.aq_find_rank <- function(tax, level) {
  # Normalize case BEFORE removing non-letters.
  norm <- function(z) {
    tolower(gsub("[^a-z]", "", tolower(as.character(z))))
  }

  level_norm <- norm(level)

  lev <- names(.aq_rank_alias)[
    vapply(
      .aq_rank_alias,
      function(a) level_norm %in% norm(a),
      logical(1)
    )
  ]

  if (length(lev) != 1L) {
    stop(
      "Unknown taxonomic level: ", level,
      ". Supported levels: ",
      paste(names(.aq_rank_alias), collapse = ", ")
    )
  }

  tax_names_norm <- norm(names(tax))
  alias_norm <- norm(.aq_rank_alias[[lev]])
  idx <- which(tax_names_norm %in% alias_norm)

  if (!length(idx)) {
    stop(
      "Taxonomy column not found for ", lev,
      ". Available columns: ",
      paste(names(tax), collapse = ", ")
    )
  }

  list(
    level = lev,
    col = names(tax)[idx[1]],
    rank = match(lev, names(.aq_rank_alias))
  )
}

#' Aggregate abundance at a taxonomic rank
#' @param data An `aq_data` object.
#' @param level Taxonomic level.
#' @param abundance `count` or `relative`.
#' @param unclassified How to handle missing labels: keep or drop.
#' @return An `aq_taxa` object.
#' @export
aq_taxa <- function(data, level = "Genus", abundance = c("count", "relative"),
                    unclassified = c("keep", "drop")) {
  if (!inherits(data, "aq_data")) stop("data must be aq_data.")
  abundance <- match.arg(abundance); unclassified <- match.arg(unclassified)
  if (is.null(data$taxonomy)) stop("Taxonomy table is required.")
  ri <- .aq_find_rank(data$taxonomy, level)
  tax <- data$taxonomy[colnames(data$otu), , drop = FALSE]
  lab <- trimws(as.character(tax[[ri$col]]))
  bad <- is.na(lab) | !nzchar(lab) | grepl("^(na|unclassified|unknown|uncultured)$", lab, ignore.case = TRUE)
  if (unclassified == "drop") {
    keep <- !bad
  } else {
    lab[bad] <- "Unclassified"
    keep <- rep(TRUE, length(lab))
  }
  if (!any(keep)) stop("No taxa remain after filtering.")
  # Use complete lineage as aggregation key to avoid merging homonymous taxa.
  parent_cols <- names(tax)[seq_len(max(0, match(ri$col, names(tax)) - 1))]
  lineage <- if (length(parent_cols)) apply(tax[, parent_cols, drop = FALSE], 1, paste, collapse = ";") else rep("", nrow(tax))
  key <- paste(lineage, lab, sep = "\034")

  # rowsum() returns one row per aggregation key and one column per sample.
  # Save those keys BEFORE transposing: after t(), rownames are sample IDs.
  summed <- rowsum(
    t(data$otu[, keep, drop = FALSE]),
    group = key[keep],
    reorder = FALSE
  )
  group_keys <- rownames(summed)
  agg <- t(summed)
  colnames(agg) <- sub(".*\034", "", group_keys)
  # Disambiguate duplicated display labels.
  if (anyDuplicated(colnames(agg))) colnames(agg) <- make.unique(colnames(agg))
  rel <- sweep(agg, 1, pmax(rowSums(agg), 1), "/")
  rownames(agg) <- rownames(data$otu); rownames(rel) <- rownames(data$otu)
  out <- structure(list(count = agg, relative = rel,
                        abundance = if (abundance == "count") agg else rel,
                        metadata = data$metadata, taxonomy = tax,
                        parameters = list(level = ri$level, abundance = abundance)),
                   class = "aq_taxa")
  out
}

#' @exportS3Method print aq_taxa
print.aq_taxa <- function(x, ...) {
  cat("AqueCo taxonomic summary:", x$parameters$level, "\n",
      nrow(x$count), "samples;", ncol(x$count), "taxa\n")
  invisible(x)
}
