\name{aq_net_plot}
\alias{aq_net_plot}
\title{Plot a microbial network with taxonomic node colors}
\description{Plots a single AqueCo network with node fill mapped to a taxonomic rank or network role. Isolate removal only changes visualization.}
\usage{
aq_net_plot(data, layout="fr", node_color="Phylum", node_size="Degree", top_n=10,
  edge_color="sign", show_legend=TRUE, remove_isolates=FALSE, seed=123)
}
\arguments{
\item{data}{Single aq_network object.}
\item{layout}{fr or circle.}
\item{node_color}{Taxonomy column, Module, Role, or none.}
\item{node_size}{Degree or Strength.}
\item{top_n}{Number of top taxonomic categories, with other categories pooled.}
\item{edge_color}{sign or none.}
\item{show_legend}{Show legends.}
\item{remove_isolates}{Hide zero-degree vertices from the plot only.}
\item{seed}{Layout random seed.}
}
