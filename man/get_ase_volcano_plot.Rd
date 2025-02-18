% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splicewiz_wrappers.R
\name{get_ase_volcano_plot}
\alias{get_ase_volcano_plot}
\title{Generate a volcano plot for diff. ASE.}
\usage{
get_ase_volcano_plot(
  x,
  sample_comparison,
  event_type = "IR",
  log2fc_cutoff = 1,
  pval_cutoff = 0.05,
  ASE_to_display = NULL,
  lab_size = 3,
  point_size = 1,
  repair_ASE = TRUE,
  col_up = "#a40000",
  col_down = "#007e2f",
  col_other = "grey",
  facet_event_type = TRUE,
  facte_nrow = NULL,
  facet_ncol = NULL,
  facet_scale = "fixed",
  facet_shrink = TRUE,
  facet_labeller = "label_value",
  facet_as.table = TRUE,
  facet_drop = TRUE,
  facet_dir = "h",
  facet_strip.position = "top",
  ...
)
}
\arguments{
\item{x}{an abject of class "parcutils_ase". This is an output of the function \code{\link[=run_ase_diff_analysis]{run_ase_diff_analysis()}}.}

\item{sample_comparison}{a character string denoting a valid differential ASE comparison. Possible comparisons can be found from x$de_comparisons.}

\item{event_type}{a character string, default "IR" denoting valid event_type to display in the volcano plot. Values could be one of the "IR","SE","AFE","ALE","MXE","A3SS","A5SS".}

\item{log2fc_cutoff}{a numeric value denoting logFC cutoff for volcano plot. Default 1.}

\item{pval_cutoff}{a numeric value denoting pvalue cutoff for volcano plot. Default 0.05.}

\item{ASE_to_display}{a character vector of the ASE to display in volcano plot, default NULL, displays non overlapping ASE.}

\item{lab_size}{a numeric value, default 3, denoting size of the labels.}

\item{point_size}{a numeric value, default 1, denoting size of the points.}

\item{repair_ASE}{logical, default TRUE, indicating whether to repair event names. See details.}

\item{col_up}{a character string, default \code{#a40000}, denoting valid color code for up regulated ASE.}

\item{col_down}{a character string, default \code{#007e2f}, denoting valid color code for up regulated ASE.}

\item{col_other}{a character string, default "grey", denoting valid color code for other than up and down regulated ASE.}

\item{facet_event_type}{logical, default TRUE, denoting whether to facet plot by event type.}

\item{facte_nrow}{pass to the argument \code{nrow} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_ncol}{pass to the argument \code{ncol} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_scale}{pass to the argument \code{scale} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_shrink}{pass to the argument \code{shrink} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_labeller}{pass to the argument \code{labeller} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_as.table}{pass to the argument \code{as.table} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_drop}{pass to the argument \code{drop} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_dir}{pass to the argument \code{dir} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{facet_strip.position}{pass to the argument \code{position} of \code{\link[ggplot2:facet_wrap]{ggplot2::facet_wrap()}}.}

\item{...}{other parameters to be passed to \code{\link[EnhancedVolcano:EnhancedVolcano]{EnhancedVolcano::EnhancedVolcano()}}.}
}
\value{
a ggplot2.
}
\description{
Generate a volcano plot for diff. ASE.
}
\details{
\itemize{
\item repair_ASE: Internally event names are taken from SpliceWiz. When repair_genes is set to TRUE the string corresponding to gene symbol will be extracted.
if one gene symbol assinged to more than one events, even names will be used instead of gene symbol (usually the case). This is useful when gene names to be revealed in the volcano plot.
}
}
\examples{
se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_pval = 1, regul_based_upon = 1)
get_ase_volcano_plot(res,  sample_comparison = res$de_comparisons[1], pval_cutoff = 1) \%>\% print()

# repair_ASE = FALSE
get_ase_volcano_plot(res,  sample_comparison = res$de_comparisons[1], pval_cutoff = 1,repair_ASE = FALSE) \%>\% print()

# ASE_to_display = "", don't show any label.
get_ase_volcano_plot(res,  sample_comparison = res$de_comparisons[1], pval_cutoff = 1,ASE_to_display = "") \%>\% print()
}
