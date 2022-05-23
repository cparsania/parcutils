# IR analysis functions
retained_introns_files  <-
  list.files(path  = "000_data_files/ir_finder_s_output/",
             pattern =  "*nondir_SRR*", full.names = T,recursive = T)

names(retained_introns_files) <- stringr::str_extract(retained_introns_files,
                                                      pattern = "SRR([^_])+")

# read irfinders output

#' Read IR-finderS output
#'
#' @param files a character vector denoting irfinderS output file(s) ending with suffix "IR-nondir".
#' @param add_prefix_chr logical, whether to add prefix 'chr' in the column seq names
#' @param remove_prefix_chr logical, whether to remove prefix 'chr' from the column seq names
#'  + chr
#'  + start
#'  + end
#'  + name
#'  + null
#'  + strand
#'  + coverage
#'  + introndepth
#'  + spliceleft
#'  + spliceright
#'  + spliceexact
#'  + irratio
#'  + warnings
#' @return
#' @export
#'
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
read_irfinderS_output <- function(files,
                                  add_prefix_chr = TRUE,
                                  remove_prefix_chr = FALSE){

  # check if each values in files are named, not null or not NA.
  stopifnot("files must be named vector" = !(is.null(files) %>% all()) | !(is.na(files) %>% any()))

  stopifnot("'add_prefix_chr' must be logical" = is.logical(add_prefix_chr))

  stopifnot("'remove_prefix_chr' must be logical" = is.logical(remove_prefix_chr))


  retained_introns_list <- purrr::map(files, ~ ..1 %>%

                                        # read
                                        readr::read_delim(delim = "\t",show_col_types = FALSE) %>%

                                        # rename to lower case
                                        dplyr::rename_all(~tolower(.)))

  # add prefix chr
  if(add_prefix_chr){
    retained_introns_list <- retained_introns_list %>%
      purrr::map( ~ ..1 %>% dplyr::mutate(chr = stringr::str_c("chr", chr, sep ="")))
  }


  # remove prefix chr
  if(remove_prefix_chr){
    retained_introns_list <- retained_introns_list %>%
      purrr::map( ~ ..1 %>% dplyr::mutate(chr = stringr::str_replace(chr, pattern = "chr",replacement = "")))
  }


  # add class
  retained_introns_list <- .assign_class_irFinderSdata(retained_introns_list)

  # add unique intron id
  retained_introns_list <- .parcutils_assign_intron_identifier(retained_introns_list)

  return(retained_introns_list)
}


#' Subset columns from IRfinder-S output
#'
#' @param x an object of class irFinderSdata
#' @param keep_columns a character vector denoting column names to keep in the output dataframe
#'  + chr
#'  + start
#'  + end
#'  + name
#'  + null
#'  + strand
#'  + coverage
#'  + introndepth
#'  + spliceleft
#'  + spliceright
#'  + spliceexact
#'  + irratio
#'  + warnings
#' @return
#' @export
#'
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' select_cols_irfinderS_output(x)
select_cols_irfinderS_output <- function(x, keep_columns = c("chr",
                                                             "start",
                                                             "end",
                                                             "name",
                                                             "null",
                                                             "strand",
                                                             "coverage",
                                                             "introndepth",
                                                             "spliceleft",
                                                             "spliceright",
                                                             "spliceexact",
                                                             "irratio",
                                                             "warnings")){

  .validate_irfinders_object(x)

  x <- purrr::map(x , ~..1 %>% dplyr::select(keep_columns))

  class(x) <- c("irFinderSdata", class(x))
  return(x)


}


#' Prepare a master list of introns
#'
#' @param f a character string denoting irfinderS output file ending with suffix "IR-nondir".
#' @param add_meta_data logical, whether to map meta data or not
#' @param bs_genome_object an object of class BSgenome
#' @param add_prefix_chr logical, whether to add prefix 'chr' in the column seqnames
#' @param remove_prefix_chr logical, whether to remove prefix 'chr' from the column seqnames
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' get_intron_master_list(f = example_files[1],add_prefix_chr = TRUE,bs_genome_object = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
get_intron_master_list <- function(f,
                                   add_meta_data = T,
                                   add_prefix_chr = FALSE,
                                   remove_prefix_chr = FALSE,
                                   bs_genome_object = NULL){


  stopifnot("'add_meta_data' must be logical." = is.logical(add_meta_data))

  names(f) <- stringr::str_c("file", length(f), sep = "_")
  x <- read_irfinderS_output(files = f,
                             add_prefix_chr = add_prefix_chr,
                             remove_prefix_chr = remove_prefix_chr)

  # if x has more than one elements, select first to prepare intron list.
  # As all elements have same rows taking any one to prepare intron list is ok.

  x <- x[[1]]

  x <- x  %>%

    # select cols
    dplyr::select(chr, start, end, name, null, strand) %>%

    # prepare gene symbol and gene id
    dplyr::mutate(gene_name = stringr::str_replace(name,"/.*" ,""),
                  gene_id = stringr::str_match(string = name,pattern = "/(.*)/")[,2],
                  intron_type = stringr::str_replace(name,".*/" ,"") ) %>%

    # prepare intron id
    dplyr::mutate(intron_id = stringr::str_c("intron_",dplyr::row_number()))

  if(add_meta_data){
    stopifnot("if add_meta_data is TRUE, bs_genome_object must be an object of class BSgenome" = is(bs_genome_object, "BSgenome"))
    colnames(x)[1] <- "seqnames"
    x <- x %>% plyranges::as_granges()
    x <- .map_granges_metadata(x = x, bs_genome_object = bs_genome_object)
    x <- x %>% tibble::as_tibble()
    colnames(x)[1] <- "chr"

  }

  return(x)
}


# add intron metadata to irfinder-s output

#' Map metadata to irfinder-s output
#' @description This function allows mapping DNA sequence, GC content and intron length to irfinder-s output
#'
#' @param x an object of class irFinderSdata
#' @param bs_genome_object an object of class BSgenome
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files[1],  add_prefix_chr = TRUE)
#' map_intron_meta_data(x = x,  bs_genome_object = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
map_intron_meta_data <- function(x , bs_genome_object = BSgenome.Hsapiens.UCSC.hg38){
  .validate_irfinders_object(x)

  stopifnot("bs_genome_object must be an object of class BSgenome" = is(bs_genome_object, "BSgenome"))

  #change name 'chr' to seqnames
  x <- purrr::map(x , function(x){
    colnames(x)[1] <- "seqnames"

    return(x)
  })

  # convert granges
  x <- purrr::map(x, ~ ..1 %>% plyranges::as_granges() )

  x_mapped <- purrr::map(x, ~ ..1 %>% .map_granges_metadata(bs_genome_object = bs_genome_object))

  x_mapped <- x_mapped %>% purrr::map(~..1 %>% tibble::as_tibble())

  # change column 'chr' to seqnames

  x_mapped <- purrr::map(x_mapped , function(x){
    colnames(x)[1] <- "chr"
    return(x)
  })
  x <- .assign_class_irFinderSdata(x_mapped)
  return(x)
}


#' Filter IR results
#' @description Based on filters applied, this function will checks whether each intron is retained or not.
#'
#' @param x an object of class irFinderSdata.
#' @param min_intron_cov a numeric value between 0 and 1, default 0.95, denoting minimum value for coverage cutoff
#' @param min_intron_depth a numeric value, default 5, denoting average sequence depth for each intron.
#' @param minimum_splice_exact a numeric value, default 5, denoting number of reads supporting intron splicing.
#' @param min_irratio a numeric value, default 0.0001, denoting ir-ratio.
#'
#' @return an object of class irFinderSdata
#' @export
#' @keywords internal
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = TRUE)
#' mark_ir_status_by_filters(x)
#'
.mark_ir_status_by_filters <- function(x ,
                                      min_intron_cov = 0.95,
                                      min_intron_depth = 5,
                                      minimum_splice_exact = 5,
                                      min_irratio = 0){
  .validate_irfinders_object(x)


  apply_intron_coverage_cutoff = TRUE

  # filter by coverage
  if(apply_intron_coverage_cutoff) {
    x_filt <- x %>%
      purrr::map(~ ..1 %>% dplyr::filter(coverage >= min_intron_cov ))
  }


  # filter by depth
  apply_depth_cutoff <- TRUE

  if(apply_depth_cutoff){
    x_filt <- x_filt %>%
      purrr::map(~ ..1 %>% dplyr::filter(introndepth >= min_intron_depth ))
  }


  # filter by splice exact
  filter_by_splice_exact <- TRUE

  if(filter_by_splice_exact){
    x_filt <- x_filt %>%
      purrr::map(~ ..1 %>% dplyr::filter(spliceexact >= minimum_splice_exact))
  }

  # filter by IR ratio
  apply_ir_ratio_cuoff <- TRUE
  if(apply_ir_ratio_cuoff){
    x_filt <- x_filt %>%
      purrr::map(~ ..1 %>%
                   dplyr::filter(irratio >= min_irratio ))
  }

  introns_remained <- x_filt %>% purrr::map(~..1 %>% dplyr::pull(intron_id))

  # map introns_remained to original object

  x <- purrr::map(names(introns_remained) , ~x[[..1]] %>%
                    dplyr::mutate(is_retained_by_filters = intron_id %in% introns_remained[[..1]]))

  names(x) <- names(x_filt)
  x <- .assign_class_irFinderSdata(x)

  return(x)

}















#' Prepare a intron reads count matrix.
#' @description Intron counts will be obtained from the column 'introndepth' of the IRFinder-S output.
#' @param x an object of class irFinderSdata
#'
#' @return a tbl.
#' @export
#'
#' @examples
#'
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' get_ir_counts(x)
#'
#'
get_ir_counts <- function(x){

  .validate_irfinders_object(x)

  y <- x %>%
    purrr::map(~ ..1 %>% dplyr::select(intron_id, "introndepth")) %>%
    tibble::enframe() %>% tidyr::unnest(cols = c(value)) %>%
    tidyr::pivot_wider(names_from = "name", values_from = "introndepth")

  return(y)

}





#' Run DEseq2 for retained introns.
#' @description This function performs DEseq2 analysis for retained introns which are identified
#' by IRFinder-S program. It mainly requires two inputs - 1) output file(s) from IRFinder-S program,
#' and 2) sample group information to inform DEseq2 which samples to be grouped for differential expression analysis.
#' Other arguments are to filter retained and differently expressed introns.
#' @param x an object of the class IRFinder-S. Usually this can be generated through the function [parcutils::read_irfinderS_output()].
#' @param min_intron_cov a value between 0 and 1, default 0.7, denoting minimum threshold for the intron coverage.
#' Introns with coverage less than the value will not be considered for DE analysis.
#' @param min_intron_depth an integer, default 10, denoting threshold for minimum number of aligned reads within an intron.
#' Same value will be passed to argument \code{min_counts} in [parcutils::run_deseq_analysis()] for differential analysis.
#' @param minimum_splice_exact an integer, default 5, denoting threshold for minimum number of spliced reads within an intron.
#' @param minimum_replicates_with_ri an integer, default 2, denoting threshold for number of replicates have an intron identified as "retained intron".
#' Same value will be passed to argument \code{min_replicates} in [parcutils::run_deseq_analysis()] for differential analysis.
#' @param sample_info internally passed to [parcutils::run_deseq_analysis()].
#' @param group_numerator internally passed to [parcutils::run_deseq_analysis()].
#' @param group_denominator internally passed to [parcutils::run_deseq_analysis()].
#' @param cutoff_lfc internally passed to [parcutils::run_deseq_analysis()].
#' @param cutoff_pval internally passed to [parcutils::run_deseq_analysis()].
#' @param cutoff_padj internally passed to [parcutils::run_deseq_analysis()].
#' @param regul_based_upon internally passed to [parcutils::run_deseq_analysis()].
#'
#' @return
#' @export
#'
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' ir_sample_info <- names(example_files) %>%
#' tibble::tibble(samples = . , groups = rep(c("control", "treatment"), each = 2) )
#' run_deseq_analysis_ir(x = x, sample_info = ir_sample_info,
#' min_intron_depth = 2, regul_based_upon =1, cutoff_pval = 0.01)
run_deseq_analysis_ir <- function(x,
                                  min_intron_cov = 0.7,
                                  min_intron_depth = 10,
                                  minimum_splice_exact = 5,
                                  minimum_replicates_with_ri = 2,
                                  sample_info,
                                  group_numerator = "treatment",
                                  group_denominator = "control",
                                  cutoff_lfc = 1,
                                  cutoff_pval = 0.05,
                                  cutoff_padj = 0.01,
                                  regul_based_upon = 1
){
  # x = x
  # min_intron_cov = 0.7
  # min_intron_depth = 2
  # minimum_splice_exact = 5
  # minimum_replicates_with_ri = 2
  # sample_info =  ir_sample_info
  # group_numerator = "treatment"
  # group_denominator = "control"
  # min_counts = 10
  # min_replicates = 2
  # cutoff_lfc = 1
  # cutoff_pval = 0.05
  # cutoff_padj = 0.01
  # regul_based_upon = 1

  .validate_irfinders_object(x)

  # mark retained introns based on filters
  x_ir_marked <- .mark_ir_status_by_filters(x,
                                           min_intron_cov = min_intron_cov,
                                           min_intron_depth = min_intron_depth,
                                           minimum_splice_exact = minimum_splice_exact)

  # prepare list of retained introns
  retained_introns_id <- x_ir_marked %>%
    purrr::map(~..1 %>% dplyr::filter(is_retained_by_filters) %>%
                 dplyr::pull(intron_id)) %>%
    tibble::enframe(name = "samples", value ="intron_id") %>%
    tidyr::unnest(cols = c(intron_id)) %>%
    dplyr::left_join(ir_sample_info, by = "samples") %>%
    dplyr::group_by(groups, intron_id) %>%
    dplyr::summarise( n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n >= minimum_replicates_with_ri) %>%
    dplyr::pull(intron_id) %>% unique()

  if(length(retained_introns_id) == 0 ){
    stop("Cannot process furhter. Number of retained introns are 0.")
  }

  # prepare intron count data for de analysis
  count_data <- get_ir_counts(x_ir_marked) %>%
    # round each value. Required to run DESeq2 analysis.
    dplyr::mutate_if(is.numeric, round,0) %>%
    # keep only those are marked as retained
    dplyr::filter(intron_id %in% retained_introns_id )

  # perform deseq analysis

  id_column <- count_data[1] %>% colnames()

  suppressWarnings(
    de_output <- run_deseq_analysis(counts = count_data,
                                    column_geneid = id_column,
                                    sample_info = sample_info,
                                    group_numerator = group_numerator,
                                    group_denominator = group_denominator,
                                    min_counts = min_intron_depth,
                                    min_replicates = minimum_replicates_with_ri,
                                    cutoff_pval = cutoff_pval,
                                    cutoff_padj = cutoff_padj ,
                                    regul_based_upon = regul_based_upon)
  )

  # prepare object parcutils_ir
  de_output <- .prepare_parcutils_ir(de_output)

  # add intron annotations to object de_output
  intron_annot <- .get_intron_annotations(x = x, intron_ids = retained_introns_id)

  # add to list column of intron_annot to final output
  de_output <- de_output %>%
    dplyr::mutate(intron_annot = list(intron_annot))

  # to have consistency across columns assign names to each elem of intron_annot
  names(de_output$intron_annot) <- de_output$comp

  return(de_output)

}


#' Annotate introns.
#'
#' @param x an object of class 'parcutils_ir'.
#' @param query_introns a character vector of intron_ids to be queried.
#' @param add_meta_data logical, default FALSE, whether to add intron meta-data (seq, GC and lenght) to the output.
#' @param bs_genome_object an object of class bs_genome. This will be used to retrieve meta-data.
#' Default  NULL. Required when add_meta_data = TRUE.
#' @param add_prefix_chr logical, default TRUE, denoting whether to add prefix 'chr' in the column seqnames.
#' @param remove_prefix_chr logical, default FALSE, denoting whether to remove prefix 'chr' from the column seqnames.
#'
#' @return a data frame
#' @export
#'
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' ir_sample_info <- names(example_files) %>%
#' tibble::tibble(samples = . , groups = rep(c("control", "treatment"), each = 2) )
#' y = run_deseq_analysis_ir(x = x, sample_info = ir_sample_info,
#' min_intron_depth = 2, regul_based_upon =1, cutoff_pval = 0.01)
#'
#' q_intron_id = c("intron_223342","intron_223728","intron_225394","intron_228167","intron_228226")
#' annotate_retained_introns(y, query_introns = q_intron_id)
annotate_retained_introns <- function(x,
                                      query_introns,
                                      add_meta_data = FALSE,
                                      bs_genome_object = NULL,
                                      add_prefix_chr = TRUE,
                                      remove_prefix_chr = FALSE
){

  # x = oo
  # query_introns = parcutils::get_genes_by_regulation(x = oo, sample_comparisons = "treatment_VS_control")
  # add_meta_data = TRUE
  # bs_genome_object = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  # add_prefix_chr = TRUE
  # remove_prefix_chr = FALSE

  .validate_parcutils_ir_object(x)

  stopifnot("'query_introns' must be a character vector." = is.character(query_introns))

  stopifnot("'add_meta_data' must be a logical." = is.logical(add_meta_data))

  stopifnot("'add_prefix_chr' must be logical." = is.logical(add_prefix_chr))

  stopifnot("'remove_prefix_chr' must be logical." = is.logical(remove_prefix_chr))

  stopifnot("'add_prefix_chr' & 'remove_prefix_chr' both cannnot be 'TRUE'." = !all(c(add_prefix_chr,remove_prefix_chr)))


  # all values in the list column 'intron_annot' carry same data. Using 1st element below.

  intron_annotations <- x$intron_annot[[1]]


  # get intron meta-data seq, GC and length
  if(add_meta_data){

    intron_annotations_mod <- intron_annotations %>%
      dplyr::rename(seqnames = chr)


    if(add_prefix_chr){
      intron_annotations_mod <- intron_annotations_mod %>%
        dplyr::mutate(seqnames = stringr::str_c("chr",seqnames, sep = ""))
    }
    if(remove_prefix_chr){
      intron_annotations_mod <- intron_annotations_mod %>%
        dplyr::mutate(seqnames = stringr::str_replace(seqnames,"chr", ""))
    }

    # convert gr
    intron_annotations_gr <- intron_annotations_mod %>%
      plyranges::as_granges()

    stopifnot("'bs_genome_object' cannot be NULL when 'add_meta_data' is TRUE")
    ir_meta_data <- .map_granges_metadata(x = intron_annotations_gr,
                                          bs_genome_object = bs_genome_object)

    # join meta-data

    intron_annotations <- intron_annotations %>%
      dplyr::left_join(ir_meta_data  %>% tibble::as_tibble() %>%
                         dplyr::select(intron_id, seq, GC,length), by = "intron_id")
  }

  # check how many queried retained introns (intron_id) not present in the all_retained

  all_retained <- intron_annotations %>%
    dplyr::pull(intron_id)

  present <- query_introns[query_introns %in% all_retained]
  not_present <- query_introns[!query_introns %in% all_retained]

  if(length(not_present) >0){
    cli::cli_warn(message =  "Intron{?s} {not_present} not found.")
  }

  if(length(present) == 0 ){
    cli::cli_abort(message =  "None of the queried introns found in the x.")
  }

  # filter and fix row order.
  out <- intron_annotations %>%
    dplyr::filter(.data$intron_id %in% present) %>%
    dplyr::slice(match(present , .data$intron_id))

  return(out)
}




#' Save differential expression results in an excel file.
#' @description All DE comparison from the \code{x} will be saved in an excel file.
#' Each tab within the excel file denotes an individual DE comparison.
#' @param x an object of class parcutils.
#' @param file a character string denoting a file name WITHOUT file extension. (e.g. "out_file" or "path/to/save/data/out_file")
#'
#' @return an absolute file path denoting a file in which data saved.
#' @export
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' ir_sample_info <- names(example_files) %>%
#' tibble::tibble(samples = . , groups = rep(c("control", "treatment"), each = 2) )
#' out <- run_deseq_analysis_ir(x = x, sample_info = ir_sample_info,
#' min_intron_depth = 2, regul_based_upon =1, cutoff_pval = 0.01)
#' save_de_results(out, file = fs::file_temp())
save_de_results <- function(x, file){

  file <- fs::path(file,ext = "xlsx")
  file <- fs::path_abs(file)
  file <- fs::file_create(path = file)

  .validate_parcutils_obj(x)

  if(inherits(x = x, what = "parcutils_ir")){
    .save_de_ir_results(x, file)
  } else{
    .save_deg_results(x, file)
  }


}







