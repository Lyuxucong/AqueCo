\name{aq_lefse_export}
\alias{aq_lefse_export}
\title{Export LEfSe differential features}
\usage{
aq_lefse_export(data, file = "LEfSe_results.csv")
}
\arguments{
\item{data}{An \code{aq_lefse} result from \code{aq_lefse()}.}
\item{file}{Output CSV file path.}
}
\value{
Invisibly returns the output file path.
}
\description{
Writes the differential-feature table stored in the
\code{res_diff} component of an \code{aq_lefse} result.
}
\examples{
\dontrun{
aq_lefse_export(lefse_result, file = "LEfSe_results.csv")
}
}
