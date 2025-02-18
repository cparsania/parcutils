% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{.perform_kemans_cluster}
\alias{.perform_kemans_cluster}
\title{Perform k-means clustering}
\usage{
.perform_kemans_cluster(data, km)
}
\arguments{
\item{data}{a dataframe. First column will be used as row identifier.
Therefore, values in the first column must be unique.}

\item{km}{an integer denoting number of clusters to return.}
}
\value{
a dataframe with two columns -1) 1st column from \code{data} and 2) clusters
}
\description{
Perform k-means clustering
}
\keyword{internal}
