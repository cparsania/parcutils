#' Find signal intensity of intergenic regions.
#' @description given the .gff file and .bw file function finds the intergenic regions and associate signal inensity of each region
#' @param bw_file a string denoting path of .bw file
#' @param gff_file a string denoting path of gff file.
#'
#' @return an object of class \link{GenomicRanges} which can be exported to various file formats such as bw, bed, begraph etc.
#' @export
#'
#' @examples
#' library("magrittr")
#' bw_file <- system.file("extdata" ,
#' "example.bw" ,
#' package = "parcutils")
#' gff_file <- system.file("extdata" ,
#'                        "C_glabrata_CBS138_version_s02-m07-r06_features.gff" ,
#' package = "parcutils")
#' oo <- parcutils::get_intergenic_singals(bw_file = bw_file,gff_file  = gff_file)
#' oo %>%
#'     rtracklayer::export.bed(con = file(paste("intergenic", ".bed",sep = "")))
#' oo %>%
#'    rtracklayer::export.bedGraph(con = file(paste("intergenic", ".bdg",sep = "")))
#' @importFrom rtracklayer import
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter mutate
#' @importFrom GenomicRanges GRanges mcolAsRleList binnedAverage gaps
#' @importFrom tools file_path_sans_ext
#' @importFrom magrittr %>%
get_intergenic_singals <- function(bw_file , gff_file){

  gff_dat <- rtracklayer::import(gff_file)

  gff_genes <- gff_dat %>%
    tibble::as_tibble() %>%
    dplyr::filter(.data$type == "gene") %>%
    GenomicRanges::GRanges()

  gff_intergenic <- gff_genes %>%
    tibble::as_tibble() %>%
    dplyr::mutate(strand = "+") %>%
    GenomicRanges::GRanges() %>%
    GenomicRanges::gaps()

  bw_dat <- rtracklayer::import(bw_file)

  bw_rtl <- GenomicRanges::mcolAsRleList(bw_dat,varname = "score")

  intergenic_ba <- GenomicRanges::binnedAverage(gff_intergenic  , bw_rtl , varname = "score")
  return(intergenic_ba)

}








