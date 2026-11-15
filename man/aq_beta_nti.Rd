\name{aq_beta_nti}
\alias{aq_beta_nti}
\title{Compute betaNTI with phylogenetic label permutations}
\description{Compute betaNTI with phylogenetic label permutations. See inst/NULLMODEL_README.md and inst/NULLMODEL_EXAMPLE.R for workflow and limitations.}
\usage{
aq_beta_nti(data, runs=500, abundance.weighted=TRUE, null.model="taxa.labels",
  seed=data$seed, workers=1)
}
\arguments{
\item{data}{Input data or grouping parameter.}
\item{runs}{See function documentation and examples.}
\item{abundance.weighted}{See function documentation and examples.}
\item{null.model}{See function documentation and examples.}
\item{seed}{See function documentation and examples.}
\item{workers}{Number of parallel workers on Unix/macOS; Windows falls back to sequential.}}
