#' Find signal intensity of intergenic regions.
#' @description given the .gff file and .bw file function finds the intergenic regions and associate signal inensity of each region
#'
#' @param bw_file a string denoting path of .bw file
#' @param gff_file a string denoting path of gff file.
#'
#' @return an object of class [GenomicRanges::granges()] which can be exported to various file formats such as bw, bed, begraph etc.
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


#' @title Convert raw read counts in to normalised values.
#' @description Reads Per Kilo Base Per Million (RPKM) and Tags Per Million (TPM) are two most waidely used normalisation methods in RNA-seq experiments.
#' This function convert raw read counts either into RPKM or TPM values.
#' @details implemntaion of RPKM and TPM can be seen in functions [parcutils::get_tpm()] and [parcutils::get_rpkm()] respectively.
#'
#' @param x A dataframe of raw counts along with mandatory columns which are \code{GeneID, Chr, Start, End, Strand, Length}.
#' @param .vars A character vector containing columns from \code{x}. Normalization will be performed only on these columns. If NULL (default) all columns normlaisation will be
#' performed on all columns.
#' @param method A character string, default TPM. Choices are one of TPM, RPKM.
#'
#' @return A dataframe with all mandatory columns along with columns mentioned in \code{.vars}. Remaining columns from \code{x} will be dropped.
#' @export
#' @import dplyr tidyselect
#' @seealso  [parcutils::get_tpm()]  [parcutils::get_rpkm()]
#' @examples
#' set.seed(123)
#' tt <- tibble::tibble(gene_id = c(paste("Gene_",1:5,sep = "")),
#'                      chr = "Chr1",
#'                      start = sample(1:100, 5),
#'                      end = sample(100:200,5),
#'                      strand = sample(c("+" ,"-"), 5, replace = TRUE),
#'                      length = (end - start )+ 1 )
#'
#' tt %<>% dplyr::mutate(sample_1 = sample(c(1:100),5)*10 ,
#'                       sample_2 = sample(c(1:100),5)*10 ,
#'                       sample_3 = sample(c(1:100),5)*100,
#'                       sample_4 = sample(c(1:100),5)*100 )
#'
#' normalise_counts(x = tt ,method = "RPKM")
#' normalise_counts(x = tt, .vars = c("sample_1","sample_2") ,method = "TPM")
#' normalise_counts(x = tt, .vars = c("sample_1","sample_2") ,method = "RPKM")
normalise_counts <- function(x, .vars = NULL ,method = "TPM"){

  x_colnames <- colnames(x)

  # validate arguments
  # x must be an object of class data.frame
  base::stopifnot(is.data.frame(x))

  # possible values for argument method : TPM, RPKM.
  base::match.arg(arg = method,choices = c("TPM","RPKM"))

  # .vars must be a character vector
  if(!(is.null(.vars) | is.character(.vars))){
    cli::cli_text("{.arg .vars} must be a {.cls NULL} or a {.cls character} vector.")
  }

  # check presence of mandatory columns in x
  mandatory_cols <- c("Gene_ID", "Chr", "Start", "End", "Strand", "Length") %>% tolower()

  if(!all(mandatory_cols %in% x_colnames)){
    i <- which(mandatory_cols %in% x_colnames == FALSE)
    cols_not_present <- mandatory_cols %>% .[i] %>% stringr::str_flatten(collapse = ",")
    cli::cli_abort("Column{?s}  {.emph {cli::col_red({cols_not_present})}} {?is/are} not present in the x.")
  }

  # select cols other than mandatory colums
  if(is.null(.vars)){
    .vars = x %>% dplyr::select(-!!mandatory_cols) %>% colnames()
  }
  .vars_quot = rlang::enquo(.vars)
  # check presence of .vars in x
  if(!all(.vars %in% x_colnames)){
    i <- which(.vars %in% x_colnames == FALSE)
    cols_not_present <- .vars %>% .[i] %>% stringr::str_flatten(collapse = ",")
    cli::cli_abort("Column{?s}  {.emph {cli::col_red({cols_not_present})}} {?is/are} not present in the x.")
  }

  # do normalization
  if(method == "TPM"){
    x %<>% dplyr::mutate_at(.vars = dplyr::vars(!!.vars_quot) ,.funs = get_tpm ,length = .$length)
  } else if(method == "RPKM") {
    x %<>% dplyr::mutate_at(.vars = dplyr::vars(!!.vars_quot) ,.funs = get_rpkm ,length = .$length)
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



#' Get alignment summary from the output of STAR program.
#' @description  STAR is widely used alignment software program for RNA-seq and other high throughput genomics data.
#' Once the alignment finishes obvious question is to know what percentage of the reads mapped to the reference.
#' These information is logged in *Log.final.out file which is one of the several files STAR generates at the end of alignment.
#' This function parse *Log.final.out file and extract the mapping statistics in a dataframe format.
#'
#' @param log_file a string denoting a path to *Log.final.out file generated by STAR aligner.
#'
#' @return a dataframe.
#' @export
#' @examples
#'
#' star_align_log_file <- system.file("extdata" , "a_Log.final.out" , package = "parcutils")
#' x =  get_star_align_log_summary(log_file = star_align_log_file)
#' print(x)
#'
#'
get_star_align_log_summary <- function(log_file){


  # check if extension of the file is valid

  file_ext <- glue::glue("{fs::path_ext_remove(log_file[[1]]) %>% fs::path_ext()}.{fs::path_ext(log_file[[1]])}")

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





#' Get intersects of upset plot.
#'
#' @param upset_data a list
#' @param upset_plot an output of [UpSetR::upset()]
#'
#' @return a tbl
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- list(A = sample(1:5, 10, replace = TRUE) , B = sample(1:5, 10, replace = TRUE) , c = sample(1:10, 10, replace = TRUE))
#' us <- UpSetR::upset(UpSetR::fromList(x))
#' get_upset_intersects(x , us )
#'
#' }
#'
get_upset_intersects <- function(upset_data, upset_plot){



  if(!inherits(x = upset_plot, "upset")){
    stop("Argument 'upset_plot' must be an object of class upset.")
  }

  if(! is(upset_data ,"list")){
    stop("Argument upset_data must be a list.")
  }


  upset_data <- unlist(upset_data, use.names = FALSE)
  upset_data <- upset_data[ !duplicated(upset_data) ]


  out <- upset_plot$New_data %>%
    tibble::as_tibble() %>%
    dplyr::mutate(us_elem = upset_data) %>%
    tidyr::gather(samples, is_element_present, -us_elem) %>%
    dplyr::filter(is_element_present == 1) %>%
    dplyr::group_by(us_elem) %>%
    tidyr::nest() %>%
    dplyr::mutate(set = purrr::map_chr(data, ~ ..1 %>% dplyr::pull(1) %>%
                                         stringr::str_c(collapse = ","))) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(set) %>%
    tidyr::nest() %>%
    dplyr::mutate(elements = purrr::map(data , ~..1 %>% dplyr::pull(1))) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup()

  return(out)

}



#' Split grouped data with names. A wrapper around dplyr::group_by()
#' ref: https://github.com/tidyverse/dplyr/issues/4223
#'
#' @param .tbl a data frame
#' @param ... arguments pass to dplyr::group_by()
#' @param keep_order logical, default TRUE, whether to maintain original order of the groups.
#'
#' @return a named list
#' @export
#'
#' @examples
#' \dontrun{
#' a <- tibble::tibble(x = 1:5, y = sample(letters[1:5]))
#' a  %>% named_group_split(y)
#' a  %>% named_group_split(y , keep_order = FALSE)
#' }
#'
#'
named_group_split <- function(.tbl, ..., keep_order = TRUE) {

  grouped <- dplyr::group_by(.tbl, ...)
  names <- rlang::inject(paste(!!!dplyr::group_keys(grouped), sep = " / "))

  gr_splt <- grouped %>%
    dplyr::group_split() %>%
    rlang::set_names(names)

  ## maintain original order of groups
  if(keep_order){
    grp_index <-  grouped %>% dplyr::group_indices() %>% unique()
    gr_splt <- gr_splt[grp_index]
  }

  return(gr_splt)
}







