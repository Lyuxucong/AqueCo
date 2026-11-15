\name{aq_null}
\alias{aq_null}
\title{Prepare filtered AqueCo null-model data}
\description{Prepare filtered AqueCo null-model data. See inst/NULLMODEL_README.md and inst/NULLMODEL_EXAMPLE.R for workflow and limitations.}
\usage{
aq_null(data, filter_thres=0.0005, group="Group", seed=123, tree=NULL, cache=FALSE)
}
\arguments{
\item{data}{Input data or grouping parameter.}
\item{filter_thres}{See function documentation and examples.}
\item{group}{Input data or grouping parameter.}
\item{seed}{See function documentation and examples.}
\item{tree}{See function documentation and examples.}
\item{cache}{Cache the pruned phylogenetic distance matrix in memory.}}
