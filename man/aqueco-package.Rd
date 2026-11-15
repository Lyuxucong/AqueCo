\name{aqueco-package}
\alias{aqueco}
\alias{aqueco-package}
\docType{package}
\title{AqueCo: Microbial Community Ecology Analysis}
\description{
AqueCo provides functions for importing, validating, analysing, and
visualizing amplicon-based microbial community data.
}
\details{
Start with \code{aq_import()} and \code{aq_check()}.
Core analyses include \code{aq_alpha()}, \code{aq_beta()},
\code{aq_taxa()}, \code{aq_ncm()}, \code{aq_rda_cca()},
and \code{aq_vpa()}. Additional modules may require suggested packages.

Bundled CSV example files are located in \code{inst/extdata}.
Use \code{system.file("extdata", "metadata.csv", package = "aqueco")}
to locate them; they are not datasets loaded by \code{data("metadata")}.
}
\seealso{
\code{\link{aq_import}}, \code{\link{aq_check}},
\code{\link{aq_alpha}}, \code{\link{aq_beta}},
\code{\link{aq_ncm}}, \code{\link{aq_rda_cca}}
}
\examples{
p <- system.file("extdata", "metadata.csv", package = "aqueco")
if (nzchar(p)) {
  md <- utils::read.csv(p, row.names = 1, check.names = FALSE)
  print(table(md$Group))
}
}
