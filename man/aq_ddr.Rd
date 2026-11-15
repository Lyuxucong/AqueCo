\name{aq_ddr}
\alias{aq_ddr}
\title{Geographic and environmental distance-decay analysis}
\usage{
aq_ddr(data, distance=c("geographic","environmental"), group=NULL, groups=NULL,
  scope=c("overall","group","both"), lon="Longitude", lat="Latitude",
  env_variables=NULL, standardize=TRUE, permutations=999, seed=123)
}
\arguments{
\item{data}{aq_data object.}
\item{distance}{geographic, environmental, or both.}
\item{group}{Metadata group column.}
\item{groups}{Optional selected group names.}
\item{scope}{overall, group, or both.}
\item{lon}{Longitude column name in metadata.}
\item{lat}{Latitude column name in metadata.}
\item{env_variables}{Numeric environmental columns; by default all numeric except coordinates.}
\item{standardize}{Standardize environmental variables before Euclidean distance.}
\item{permutations}{Mantel permutation count.}
\item{seed}{Random seed.}
}
\description{Fits pairwise Bray-Curtis against geographic haversine or environmental Euclidean distances. Reports descriptive linear R-squared and Mantel permutation P.}
