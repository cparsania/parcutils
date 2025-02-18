% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splicewiz_wrappers.R
\name{get_diff_ASE_count_barplot}
\alias{get_diff_ASE_count_barplot}
\title{Generate a barplot of diff ASE counts.}
\usage{
get_diff_ASE_count_barplot(
  x,
  col_up = "#a40000",
  col_down = "#16317d",
  font_size = 12,
  event_type = NULL,
  sample_comparisons = NULL,
  show_counts = TRUE,
  ...
)
}
\arguments{
\item{x}{an object of class parcutils_ase.}

\item{col_up}{a character string, default \verb{#a40000}, denoting valid a color name for "Up" regulated genes.}

\item{col_down}{a character string, default \verb{#16317d}, denoting valid a color name "Down" regulated genes.}

\item{font_size}{a numeric, default 12, denoting font size in the plot.}

\item{event_type}{a character vector or a string denoting type of events to plot. Default NULL, includes all.}

\item{sample_comparisons}{a character vector or a string denoting sample comparisons to plot. Default NULL, includes all.}

\item{show_counts}{a logical, default \code{TRUE}, denoting whether to show counts on each bar.}

\item{...}{Other arguments pass to the function \code{\link[ggplot2:facet_grid]{ggplot2::facet_grid()}}}
}
\value{
a bar plot.
}
\description{
Generate a barplot of diff ASE counts.
}
\examples{
se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_padj = 1, regul_based_upon = 2)
SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
get_diff_ASE_count_barplot(res, event_type = c("IR") )
}
