\name{aq_lefse}
\alias{aq_lefse}
\title{LEfSe analysis powered by microeco}
\usage{
aq_lefse(data, group, min_total_abundance = 10, alpha = 0.05,
  p_adjust_method = "none", ...)
}
\arguments{
\item{data}{An aq_data object with taxonomy and metadata.}
\item{group}{Metadata grouping column.}
\item{min_total_abundance}{Minimum total count per ASV.}
\item{alpha}{Significance cutoff.}
\item{p_adjust_method}{P-value adjustment method passed to microeco.}
\item{...}{Additional microeco trans_diff arguments.}
}
\description{Converts AqueCo data to microeco::microtable and calls trans_diff with method lefse.}
