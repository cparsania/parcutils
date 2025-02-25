.onLoad <- function(libname, pkgname) {
  required_pkgs <- c("plyranges",
                     "enrichplot",
                     "EnrichedHeatmap",
                     "EnhancedVolcano",
                     "DESeq2",
                     "ComplexHeatmap",
                     "clusterProfiler",
                     "BSgenome",
                     "rtracklayer",
                     "org.Mm.eg.db",
                     "GenomicRanges",
                     "AnnotationDbi")

  missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]

  if (length(missing_pkgs) > 0) {
    message("Installing missing Bioconductor dependencies: ", paste(missing_pkgs, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(missing_pkgs)
  }
}
