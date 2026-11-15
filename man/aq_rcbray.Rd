\name{aq_rcbray}
\alias{aq_rcbray}
\title{Compute abundance-informed RCbray}
\description{Compute abundance-informed RCbray. See inst/NULLMODEL_README.md and inst/NULLMODEL_EXAMPLE.R for workflow and limitations.}
\usage{
aq_rcbray(data, runs=1000, seed=data$seed, workers=1)
}
\arguments{
\item{data}{Input data or grouping parameter.}
\item{runs}{See function documentation and examples.}
\item{seed}{See function documentation and examples.}
\item{workers}{Number of parallel workers on Unix/macOS; Windows falls back to sequential.}}
