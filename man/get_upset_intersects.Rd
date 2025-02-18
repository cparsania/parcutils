% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_bioinfo.R
\name{get_upset_intersects}
\alias{get_upset_intersects}
\title{Get intersects of upset plot.}
\usage{
get_upset_intersects(upset_data, upset_plot)
}
\arguments{
\item{upset_data}{a list}

\item{upset_plot}{an output of \code{\link[UpSetR:upset]{UpSetR::upset()}}}
}
\value{
a tbl
}
\description{
Get intersects of upset plot.
}
\examples{
\dontrun{
set.seed(123)
x <- list(A = sample(1:5, 10, replace = TRUE) , B = sample(1:5, 10, replace = TRUE) , c = sample(1:10, 10, replace = TRUE))
us <- UpSetR::upset(UpSetR::fromList(x))
get_upset_intersects(x , us )

}

}
