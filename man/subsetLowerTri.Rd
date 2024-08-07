% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/subsetLowerTriR.R
\name{subsetLowerTri}
\alias{subsetLowerTri}
\title{subsetLowerTri}
\usage{
subsetLowerTri(m)
}
\arguments{
\item{m}{distance matrix}
}
\value{
vector of elements in the lower triangle of the distance matrix (input m)
}
\description{
This function is within the code for finding correlation of distances. It takes a distance matrix as an input,
and returns the lower triangle. This is useful because a distance matrix is symmetric, so we can take the lower
triangle only to reduce computation time
}
\keyword{internal}
