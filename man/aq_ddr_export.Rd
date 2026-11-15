\name{aq_ddr_export}
\alias{aq_ddr_export}
\title{Export distance-decay analysis results}
\usage{
aq_ddr_export(data, prefix = "DDR")
}
\arguments{
\item{data}{An \code{aq_ddr} result from \code{aq_ddr()}.}
\item{prefix}{Output file prefix, optionally including a directory path.}
}
\value{
Invisibly returns the paths of two CSV files.
}
\description{
Exports the pairwise distances and the model summary from a
distance-decay analysis.
}
\details{
Writes \code{<prefix>_pairs.csv} and \code{<prefix>_summary.csv}.
The destination directory must already exist.
}
\examples{
\dontrun{
aq_ddr_export(ddr_result, prefix = "results/DDR")
}
}
