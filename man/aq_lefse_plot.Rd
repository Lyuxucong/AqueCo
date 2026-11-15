\name{aq_lefse_plot}
\alias{aq_lefse_plot}
\title{Plot LEfSe LDA bars with controlled feature counts or a cladogram}
\usage{
aq_lefse_plot(data, type=c("lda","cladogram"), threshold=2, group_order=NULL,
  colors=NULL, top_n=NULL, top_n_per_group=NULL, tax_level=NULL, short_labels=TRUE,
  ...)
}
\arguments{
\item{data}{aq_lefse result.}
\item{type}{lda or cladogram.}
\item{threshold}{Minimum absolute LDA score for plotting.}
\item{group_order}{Optional group order.}
\item{colors}{Group color vector, preferably named.}
\item{top_n}{Maximum total features displayed.}
\item{top_n_per_group}{Maximum features per group.}
\item{tax_level}{Optional rank: Kingdom, Phylum, Class, Order, Family, Genus, Species. Requires terminal rank prefixes in Taxa.}
\item{short_labels}{Show terminal rank label rather than full lineage.}
\item{...}{For lda: width and use_number (alias for top_n); for cladogram: microeco arguments.}
}
\description{Filters LEfSe result rows before plotting LDA bars. This does not rerun statistical tests.}
