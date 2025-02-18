% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_related.R
\name{gseMsigDB}
\alias{gseMsigDB}
\title{Perform gene set enrichment analysis against Molecular Signature Database (MSigDB).}
\usage{
gseMsigDB(
  gene_list,
  from_type = "SYMBOL",
  orgdb = org.Mm.eg.db::org.Mm.eg.db,
  msigdb_category = "H",
  species = "Mus musculus",
  ...
)
}
\arguments{
\item{gene_list}{a named numeric vector. Values in the vector will be used to rank the genes in GSEA analysis. Values can be normalised expression, log2FC, test statistics etc. Vector names denote gene names.}

\item{from_type}{a character string, default "SYMBOL", denotes type of gene identifier. Possible values could be one of those available as keys of OrgDB object.}

\item{orgdb}{an object of the class OrgDB. Default org.Mm.eg.db::org.Mm.eg.db}

\item{msigdb_category}{a character string, default "H", denoting one of the nine catagories of geneList given in the msigdb. Possible values are: "H","C1","C2","C3","C4","C5","C6","C7", and "C8".}

\item{species}{a character string denoting name of the species. Default "Mus musculus". Possible values can be found via the function \code{\link[msigdbr:msigdbr_species]{msigdbr::msigdbr_species()}}.}

\item{...}{other arguments to be passed to \code{\link[clusterProfiler:GSEA]{clusterProfiler::GSEA()}}.}
}
\value{
an output of \code{\link[clusterProfiler:GSEA]{clusterProfiler::GSEA()}}.
}
\description{
This is a wrapper around \code{\link[clusterProfiler:GSEA]{clusterProfiler::GSEA()}} to perform GSEA analysis against MsigDB.
}
\examples{
\dontrun{

}
}
