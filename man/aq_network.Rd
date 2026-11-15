\name{aq_network}
\alias{aq_network}
\title{AqueCo network analysis: aq_network}
\description{Network analysis function; see inst/NETWORK_README.md for details and limitations.}
\usage{
aq_network(data, group="Group", groups=NULL, method="spearman",
  min_abundance=0.0005, min_prevalence=0.3, cor_cutoff=0.6, p_cutoff=0.05,
  p_adjust="BH", min_samples=10)
}
\arguments{
\item{group}{Input data or grouping parameter.}
\item{data}{See module documentation in inst/NETWORK_README.md.}
\item{groups}{See module documentation in inst/NETWORK_README.md.}
\item{method}{See module documentation in inst/NETWORK_README.md.}
\item{min_abundance}{See module documentation in inst/NETWORK_README.md.}
\item{min_prevalence}{See module documentation in inst/NETWORK_README.md.}
\item{cor_cutoff}{See module documentation in inst/NETWORK_README.md.}
\item{p_cutoff}{See module documentation in inst/NETWORK_README.md.}
\item{p_adjust}{See module documentation in inst/NETWORK_README.md.}
\item{min_samples}{See module documentation in inst/NETWORK_README.md.}
}
