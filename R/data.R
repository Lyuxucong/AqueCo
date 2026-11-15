#' Construct an AqueCo data object
#' @param otu Numeric abundance matrix; rows are samples and columns are taxa.
#' @param taxonomy Taxonomy data.frame with taxa IDs as row names.
#' @param metadata Sample metadata with sample IDs as row names.
#' @param env Optional environmental data, with sample IDs as row names.
#' @param tree Optional phylo tree.
#' @return An `aq_data` object.
#' @export
aq_data <- function(otu, taxonomy = NULL, metadata = NULL, env = NULL, tree = NULL) {
  otu <- as.matrix(otu)
  if (!is.numeric(otu)) stop("otu must be numeric.")
  if (is.null(rownames(otu)) || is.null(colnames(otu)))
    stop("otu must have sample row names and taxa column names.")
  if (anyDuplicated(rownames(otu)) || anyDuplicated(colnames(otu)))
    stop("Sample IDs and taxa IDs must be unique.")
  if (any(!is.finite(otu)) || any(otu < 0))
    stop("otu must contain finite, non-negative values.")
  if (!is.null(taxonomy)) {
    taxonomy <- as.data.frame(taxonomy, check.names = FALSE)
    if (is.null(rownames(taxonomy))) stop("taxonomy row names must be taxa IDs.")
  }
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata, check.names = FALSE)
    if (is.null(rownames(metadata))) stop("metadata row names must be sample IDs.")
    if (!all(rownames(otu) %in% rownames(metadata)))
      stop("metadata is missing one or more OTU sample IDs.")
    metadata <- metadata[rownames(otu), , drop = FALSE]
  }
  if (!is.null(env)) {
    env <- as.data.frame(env, check.names = FALSE)
    if (is.null(rownames(env))) stop("env row names must be sample IDs.")
    if (!all(rownames(otu) %in% rownames(env)))
      stop("env is missing one or more OTU sample IDs.")
    env <- env[rownames(otu), , drop = FALSE]
  }
  structure(list(otu = otu, taxonomy = taxonomy, metadata = metadata,
                 env = env, tree = tree),
            class = "aq_data")
}

#' Import data tables into AqueCo
#'
#' The input OTU/ASV table must already have taxon IDs as row names
#' (for CSV files, use read.csv(..., row.names = 1)). Rows are taxa
#' and columns are samples. The table is transposed internally so that
#' samples are stored in rows and taxa in columns. No first-column
#' detection or removal is performed.
#'
#' @param otu OTU/ASV abundance matrix or data.frame already read with
#'   taxon IDs in row names (for CSV input, use `row.names = 1`).
#'   Rows must be taxa and columns must be samples; no ID column is detected
#'   or removed by `aq_import()`.
#' @param taxonomy Taxonomy table; taxa IDs should be row names.
#' @param metadata Sample metadata; sample IDs should be row names.
#' @param env Environmental data; sample IDs should be row names.
#' @param tree Optional phylo object.
#' @return An `aq_data` object with samples in rows and taxa in columns.
#' @export
aq_import <- function(otu, taxonomy = NULL, metadata = NULL, env = NULL,
                      tree = NULL) {

  # Input is already read with row.names = 1:
  # rows = taxa, columns = samples.
  if (is.data.frame(otu)) {
    otu <- as.matrix(otu)
  } else if (!is.matrix(otu)) {
    stop("otu must be a matrix or data.frame read with row.names = 1.")
  }

  if (is.null(rownames(otu)) || is.null(colnames(otu))) {
    stop(
      "OTU table must have taxon IDs as row names and sample IDs as column names. ",
      "Please read the OTU file with read.csv(..., row.names = 1)."
    )
  }

  if (anyNA(rownames(otu)) || any(!nzchar(trimws(rownames(otu)))) ||
      anyDuplicated(rownames(otu))) {
    stop("OTU/ASV row names must be unique and non-empty.")
  }

  if (anyNA(colnames(otu)) || any(!nzchar(trimws(colnames(otu)))) ||
      anyDuplicated(colnames(otu))) {
    stop("OTU sample column names must be unique and non-empty.")
  }

  suppressWarnings(storage.mode(otu) <- "numeric")

  if (anyNA(otu) || any(!is.finite(otu)) || any(otu < 0)) {
    stop("OTU abundance must contain finite, non-negative numeric values.")
  }

  # External format: Taxa x Samples.
  # Internal format: Samples x Taxa.
  otu <- t(otu)

  aq_data(
    otu = otu,
    taxonomy = taxonomy,
    metadata = metadata,
    env = env,
    tree = tree
  )
}

#' Check an AqueCo data object
#' @param x An `aq_data` object.
#' @return Invisibly returns a diagnostic data.frame.
#' @export
aq_check <- function(x) {
  if (!inherits(x, "aq_data")) stop("x must be an aq_data object.")
  issues <- character()
  if (any(rowSums(x$otu) == 0)) issues <- c(issues, "One or more samples have zero library size.")
  if (any(colSums(x$otu) == 0)) issues <- c(issues, "One or more taxa have zero total abundance.")
  if (!is.null(x$taxonomy)) {
    miss <- setdiff(colnames(x$otu), rownames(x$taxonomy))
    if (length(miss)) issues <- c(issues, paste(length(miss), "taxa lack taxonomy annotations."))
  }
  if (!is.null(x$metadata) && !identical(rownames(x$metadata), rownames(x$otu)))
    issues <- c(issues, "Metadata order differs from OTU sample order.")
  out <- data.frame(status = if (length(issues)) "warning" else "ok",
                    message = if (length(issues)) issues else "No basic data integrity issues detected.")
  print(out, row.names = FALSE)
  invisible(out)
}

#' @exportS3Method print aq_data
print.aq_data <- function(x, ...) {
  cat("AqueCo data object\nSamples:", nrow(x$otu), "\nTaxa:", ncol(x$otu),
      "\nTaxonomy:", !is.null(x$taxonomy), "\nMetadata:", !is.null(x$metadata), "\n")
  invisible(x)
}
