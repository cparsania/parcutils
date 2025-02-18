% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/splicewiz_wrappers.R
\name{get_ASE_by_regulation}
\alias{get_ASE_by_regulation}
\alias{get_ASEsets_by_regulation}
\title{Get Alternate Splice Events based on their differential regulation.}
\usage{
get_ASE_by_regulation(
  x,
  sample_comparisons,
  regulation = "both",
  event_type = NULL,
  simplify = FALSE
)

get_ASEsets_by_regulation(
  x,
  sample_comparisons,
  regulation = "both",
  event_type = NULL
)
}
\arguments{
\item{x}{an abject of class "parcutils_ase". This is an output of the function \code{\link[=run_ase_diff_analysis]{run_ase_diff_analysis()}}.}

\item{sample_comparisons}{a character vector denoting  sample comparisons for which ASE to be obtained.}

\item{regulation}{a character string, default \code{both}. Values can be one of the \code{up}, \code{down}, \code{both}, \code{other}, \code{all}.
\itemize{
\item \code{up} : returns all up regulated ASE.
\item \code{down} : returns all down regulated ASE
\item \code{both} : returns all up and down regulated ASE.
\item \code{other} : returns ASE other than up and down regulated ASE.
\item \code{all} : returns all ASE.
}}

\item{event_type}{a character string denoting event type for which one of the \code{up}, \code{down}, \code{both}, \code{other} or \code{all} events to be returned based.
Choices are: \code{IR}, \code{SE}, \code{AFE}, \code{ALE}, \code{A3SS}, \code{A5SS} and \code{MXE.} Default NULL, uses all event types.}

\item{simplify}{logical, default FALSE, if TRUE returns result in a dataframe format.}
}
\value{
a list or dataframe.
}
\description{
Get Alternate Splice Events based on their differential regulation.
}
\details{
For a given differential comparison this function returns \code{up}, \code{down}, \code{both}, \code{other} and \code{all} ASE.
}
\examples{
se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)

res <- run_ase_diff_analysis(x = se, test_factor = "replicate", test_nom = c("P","Q","R") ,test_denom = c("Q","R","P"),  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_pval = 1, regul_based_upon = 1)

# get both up and down regulated genes
get_ASE_by_regulation(x = res, sample_comparisons = c("P_VS_Q")) \%>\% head()

# get up genes only
get_ASE_by_regulation(x = res, sample_comparisons = c("P_VS_Q") , regul = "up") \%>\% head()

# get SE events
get_ASE_by_regulation(x = res, sample_comparisons = c("P_VS_Q") , regul = "down",event_type = "SE") \%>\% head()

# get genes other than up and down
get_ASE_by_regulation(x = res, sample_comparisons = c("P_VS_Q") , regul = "other") \%>\% head()

# Simplify output for multiple sample comparisons
get_ASE_by_regulation(x = res, sample_comparisons = res$de_comparisons, simplify = TRUE, regul= "both",event_type = "IR")


# get genesets by regulation. It uses sample comparison and regulation to name each output geneset.

get_ASEsets_by_regulation(x = res, sample_comparisons = "P_VS_Q", regul = "both")

}
