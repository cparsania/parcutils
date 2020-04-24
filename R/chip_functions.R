#' Find signal intensity of intergenic regions.
#' @description given the .gff file and .bw file function finds the intergenic regions and associate signal inensity of each region
#' @param bw_file a string denoting path of .bw file
#' @param gff_file a string denoting path of gff file.
#'
#' @return an object of class \link{GenomicRanges} which can be exported various file formats such as bw, bed, begraph etc.
#' @export
#'
#' @examples
#'
#'
#'
get_intergenic_singals <- function(bw_file , gff_file){


  gff_dat <- rtracklayer::import(gff_file)

  gff_genes <- gff_dat %>%
    tibble::as_tibble() %>%
    dplyr::filter(type == "gene") %>%
    GenomicRanges::GRanges()

  gff_intergenic <- gff_genes %>%
    tibble::as_tibble() %>%
    dplyr::mutate(strand = "+") %>%
    GenomicRanges::GRanges() %>%
    gaps()

  bw_dat <- rtracklayer::import(bw_file)

  bw_rtl <- GenomicRanges::mcolAsRleList(bw_dat,varname = "score")

  intergenic_ba <- GenomicRanges::binnedAverage(gff_intergenic  , bw_rtl , varname = "score")
  return(intergenic_ba)

}

bw_files <- list.files("/Users/chiragparsania/Documents/Projects/10_Nandan/iMac_backup/Cg_Pol_II/2_fastq_and_alignment/xbp1_deletion/180518_D00691_0103_BCCHBVANXX/mapping/" ,
                       pattern = "*.bw" , recursive = T ,full.names = T)

gff_file <- "/Users/chiragparsania/Documents/Projects/10_Nandan/iMac_backup/Cg_Pol_II/1_genome/CBS138_s02-m07-r06/annotation/C_glabrata_CBS138_version_s02-m07-r06_features.gff"

bw_file <- bw_files[1]
bw_file_name <-  tools::file_path_sans_ext(bw_file)

oo <- get_intergenic_singals(bw_file = bw_file,gff_file  = gff_file)

oo %>% export.bed(con = file(paste(bw_file_name,"_intergenic", ".bed",sep = "")))
oo %>% export.bedGraph(con = file(paste(bw_file_name,"_intergenic", ".bdg",sep = "")))







