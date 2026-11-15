#' AqueCo: Microbial Community Ecology Analysis
#'
#' AqueCo provides functions for importing, validating, analysing, and
#' visualizing amplicon-based microbial community data.
#'
#' Start with [aq_import()] and [aq_check()]. Core analyses include
#' [aq_alpha()], [aq_beta()], [aq_taxa()], [aq_ncm()],
#' [aq_rda_cca()], and [aq_vpa()]. Optional modules may require
#' packages listed in `Suggests`.
#'
#' Bundled CSV example files are located in `inst/extdata`.
#' Use `system.file("extdata", "metadata.csv", package = "aqueco")`
#' to locate them. They are not datasets loaded by `data("metadata")`.
#'
#' @name aqueco-package
#' @aliases aqueco
#' @docType package
#' @examples
#' p <- system.file("extdata", "metadata.csv", package = "aqueco")
#' if (nzchar(p)) {
#'   md <- utils::read.csv(p, row.names = 1, check.names = FALSE)
#'   print(table(md$Group))
#' }
NULL
