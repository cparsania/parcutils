% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{.generate_go_emap_plot}
\alias{.generate_go_emap_plot}
\title{Generate a list of GO EMAP plot.}
\usage{
.generate_go_emap_plot(x, n_terms = 30, color_terms_by = "p.adjust")
}
\arguments{
\item{x}{a list containing output of GO enrichment. Each element of the list is an output of  \code{\link[clusterProfiler:enrichGO]{clusterProfiler::enrichGO()}}.}

\item{n_terms}{number of terms to show.}

\item{color_terms_by}{}
}
\value{
a list of GO plots.
}
\description{
Generate a list of GO EMAP plot.
}
\keyword{internal}
