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
#' \dontrun{
#' #' library("magrittr")
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
#' }
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
#' @param x A dataframe of raw counts along with mandatory columns which are \code{GeneID, Chr, Start, End, Strand, Length}.
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
#' tt <- tibble::tibble(GeneID = c(paste("Gene_",1:5,sep = "")),
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
  mandatory_cols <- c("GeneID", "Chr", "Start", "End", "Strand", "Length")

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



#' Get alignment summary of STAR aligner
#'
#' @param log_file a string denoting a path to *Log.final.out file generated by STAR aligner.
#'
#' @return a tbl.
#' @import fs
#' @export
get_star_aling_log_summary <- function(log_file){


  # check if extension of the file is valid

  file_ext <- glue::glue("{fs::path_ext_remove(star_log_files[[1]]) %>% fs::path_ext()}.{fs::path_ext(star_log_files[[1]])}")

  if(file_ext != "final.out") {
    stop("Extension of file must be \"*Log.final.out")
  }

  # read file
  file_cont <- log_file %>% readr::read_lines() %>% stringr::str_trim()

  ## Each element of the vector is a line from *Log.final.out file. get required data from each element of the vector

  #1                                  Started job on |	Jan 11 10:06:42
  #2                             Started mapping on |	Jan 11 10:13:03
  #3                                     Finished on |	Jan 11 10:32:12
  #4        Mapping speed, Million of reads per hour |	197.65
  #5
  #6                           Number of input reads |	63084121
  #7                       Average input read length |	300
  #8                                     UNIQUE READS:
  #9                    Uniquely mapped reads number |	60272060
  #10                         Uniquely mapped reads % |	95.54%
  #11                           Average mapped length |	298.09
  #12                        Number of splices: Total |	61018038
  #13             Number of splices: Annotated (sjdb) |	60401463
  #14                        Number of splices: GT/AG |	60461428
  #15                        Number of splices: GC/AG |	421402
  #16                        Number of splices: AT/AC |	62147
  #17                Number of splices: Non-canonical |	73061
  #18                       Mismatch rate per base, % |	0.23%
  #19                          Deletion rate per base |	0.01%
  #20                         Deletion average length |	1.88
  #21                         Insertion rate per base |	0.01%
  #22                        Insertion average length |	1.64
  #23                              MULTI-MAPPING READS:
  #24         Number of reads mapped to multiple loci |	1309690
  #25              % of reads mapped to multiple loci |	2.08%
  #26         Number of reads mapped to too many loci |	13337
  #27              % of reads mapped to too many loci |	0.02%
  #28                                   UNMAPPED READS:
  #29   Number of reads unmapped: too many mismatches |	0
  #30        % of reads unmapped: too many mismatches |	0.00%
  #31             Number of reads unmapped: too short |	1468446
  #32                  % of reads unmapped: too short |	2.33%
  #33                 Number of reads unmapped: other |	20588
  #34                      % of reads unmapped: other |	0.03%
  #35                                   CHIMERIC READS:
  #36                        Number of chimeric reads |	0
  #37                             % of chimeric reads |	0.00%



  elems_to_keep <- c(total_input_reads = 6,
                     reads_mapped_uniquely = 9 ,
                     prcnt_mapped_uniquely = 10,
                     avg_input_read_length = 7 ,
                     avg_mapped_length = 11,
                     reads_mapped_to_multi_loci = 24,
                     prcnt_mapped_to_multi_loci = 25,
                     reads_mapped_to_many_loci =  26,
                     prcnt_mapped_to_many_loci =  27,
                     reads_unmapped_too_short = 31,
                     prcnt_unmapped_too_short = 32,
                     reads_unmapped = 33,
                     prcnt_unmapped = 34)

  # convert in to a table
  file_cont[elems_to_keep] %>%
    stringr::str_split(pattern = " \\|\\t") %>%
    purrr::map_df(~ tibble::tibble(type = ..1[1], val = ..1[2]))

}






