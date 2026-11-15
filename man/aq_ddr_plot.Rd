\name{aq_ddr_plot}
\alias{aq_ddr_plot}
\title{Plot distance-decay relationships}
\usage{
aq_ddr_plot(data, distance = NULL, groups = NULL, compare = TRUE,
  show_ci = TRUE, show_points = TRUE)
}
\arguments{
\item{data}{An \code{aq_ddr} result from \code{aq_ddr()}.}
\item{distance}{Optional subset: \code{"geographic"} and/or
\code{"environmental"}. \code{NULL} displays available distances.}
\item{groups}{Optional group names to display.}
\item{compare}{If \code{TRUE}, overlay group fits; otherwise facet by group.}
\item{show_ci}{If \code{TRUE}, show linear-model confidence ribbons.}
\item{show_points}{If \code{TRUE}, show sample-pair points.}
}
\value{
A \code{ggplot} object.
}
\description{
Visualizes pairwise community dissimilarity against geographic or
environmental distance using results from \code{aq_ddr()}.
}
\details{
Confidence ribbons are based on descriptive linear fits; they are not
confidence intervals for the Mantel permutation test.
}
\examples{
\dontrun{
print(aq_ddr_plot(ddr_result))
}
}
