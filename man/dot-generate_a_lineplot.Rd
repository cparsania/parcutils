% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{.generate_a_lineplot}
\alias{.generate_a_lineplot}
\title{Generate a gene expression line plot.}
\usage{
.generate_a_lineplot(
  line_plot_data_wide,
  km,
  facet_clusters,
  scale_log10,
  line_transparency,
  show_average_line,
  average_line_color,
  average_line_size,
  average_line_summary_method
)
}
\arguments{
\item{line_plot_data_wide}{a dataframe containing data for line plot.}

\item{km}{a numeric, indicates number of clusters for k-means clustering.}

\item{facet_clusters}{a logical, denotes whether to facet clusters generated by k-means.}

\item{scale_log10}{a logical, denotes whether to transform Y axis on log10 scale.}

\item{line_transparency}{a numeric, denotes value for line transparency.}

\item{show_average_line}{a logical, denotes whether to show average line.}

\item{average_line_color}{a character string, denotes color for average line.}

\item{average_line_size}{a numeric, denotes size for average line.}

\item{average_line_summary_method}{character string either \code{mean} or  \code{median},
denotes summary method for average line.}
}
\value{
a ggplot2.
}
\description{
This function called internally from the function
\code{get_gene_expression_line_plot}.
}
\keyword{internal}
