#' Find signal intensity of intergenic regions.
#' @description given the .gff file and .bw file function finds the intergenic regions and associate signal inensity of each region
#'
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
#' oo <- parcutils::get_intergenic_signals(bw_file = bw_file,gff_file  = gff_file)
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
get_intergenic_signals <- function(bw_file , gff_file){

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


#' @title Derive normalized counts from raw read numbers
#' @description This function allows converting raw read counts into normalized counts either by method Tags Per Million (TPM) or Reads Per Kilo-base Per Million (RPKM)
#' @details implemntaion of RPKM and TPM can be seen in functions \link[parcutils]{get_tpm} and \link[parcutils]{get_rpkm} respectively.
#'
#' @param x A dataframe of raw counts along with mandatory columns which are \code{Geneid, Chr, Start, End, Strand, Length}.
#' @param .vars A character vector containing columns from \code{x}. Normalization will be performed on these columns.
#' @param method A character string, default TPM. Choices are one of TPM, RPKM.
#'
#' @return A dataframe with all mandatory columns along with columns mentioned in \code{.vars}. Remaining columns from \code{x} will be dropped.
#' @export
#' @import dplyr tidyselect
#' @importFrom glue glue
#' @examples
#' \dontrun{
#' set.seed(123)
#' tt <- tibble::tibble(Geneid = c(paste("Gene_",1:5,sep = "")),
#'                      Chr = "Chr1",
#'                      Start = sample(1:100, 5),
#'                      End = sample(100:200,5),
#'                      Strand = sample(c("+" ,"-"), 5, replace = T),
#'                      Length = (End - Start )+ 1 )
#'
#' tt %<>% dplyr::mutate(sample_1 = sample(c(1:100),5)*10 ,
#'                       sample_2 = sample(c(1:100),5)*10 ,
#'                       sample_3 = sample(c(1:100),5)*100,
#'                       sample_4 = sample(c(1:100),5)*100 )
#'
#'
#' normalise_counts(x = tt, .vars = c("sample_1","sample_2") ,method = "TPM")
#' normalise_counts(x = tt, .vars = c("sample_1","sample_2") ,method = "RPKM")
#' }
normalise_counts <- function(x, .vars ,method = "TPM"){

  x_colnames <- colnames(x)
  .vars_quot = rlang::enquo(.vars)

  # validate arguments
  # x must be an object of class data.frame
  base::stopifnot(is.data.frame(x))

  # possible values for argument method : TPM, RPKM.
  base::match.arg(arg = method,choices = c("TPM","RPKM"))

  # .vars must be a character vector
  base::stopifnot(is.character(.vars))

  # check presence of mandatory columns in x
  mandatory_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

  if(!all(mandatory_cols %in% x_colnames)){
    i <- which(mandatory_cols %in% x_colnames == F)
    cols_not_present <- mandatory_cols %>% .[i] %>% stringr::str_flatten(collapse = ",")
    stop(glue::glue("Columns {cols_not_present} are not present in the x"))
  }

  # check presence of .vars in x
  if(!all(.vars %in% x_colnames)){
    i <- which(.vars %in% x_colnames == F)
    cols_not_present <- .vars %>% .[i] %>% stringr::str_flatten(collapse = ",")
    stop(glue::glue("Columns {cols_not_present} are not present in the x"))
  }

  # do normalization
  if(method == "TPM"){
    x %<>% dplyr::mutate_at(.vars = dplyr::vars(!!.vars_quot) ,.funs = get_tpm ,length = .$Length)
  } else if(method == "RPKM") {
    x %<>% dplyr::mutate_at(.vars = dplyr::vars(!!.vars_quot) ,.funs = get_rpkm ,length = .$Length)
  } else{
    stop("arg should be one  either TPM or RPKM.")
  }

  x %<>% dplyr::select(tidyselect::all_of(mandatory_cols), !!.vars_quot)

  return(x)
}


#' @title Calculate RPKM (Reads per kilobase per million).
#' @description given the raw counts and length of genomic features this function returns normalized counts as reads per kilo-base per million  (RPKM).
#' @param counts a numeric vector giving counts of reads aligned against genomic intervals.
#' @param lengths a numeric vector giving length of genomic intervals.
#'
#' @return a numeric vector giving RPKM
#' @keywords internal
#'
get_rpkm <- function(counts, lengths) {
  stopifnot(is.numeric(counts))
  stopifnot(is.numeric(lengths))

  # normalize to length
  rate <- counts / lengths
  # normalize to total counts (total mapped reads)
  rate / sum(counts) * 1e6
}


#' @title Calculate TPM (Tags Per Million).
#' @description given the raw counts and length of genomic features this function returns normalized counts as tags per million (TPM).
#' @param counts a numeric vector giving counts of reads aligned against genomic intervals.
#' @param lengths a numeric vector giving length of genomic intervals.
#'
#' @return a numeric vector giving TPM
#' @keywords internal
get_tpm <- function(counts, lengths) {
  stopifnot(is.numeric(counts))
  stopifnot(is.numeric(lengths))
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}





