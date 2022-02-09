#' @title Perform differential expression analysis using [DESeq2::DESeq()]
#' @description
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This is a wrapper function build upon [DESeq2::DESeq()]  and  [DESeq2::DESeqResults()].
#' to find diff genes and categories them based on various cutoffs such p-value, padj-value, log2fc etc.
#' It also allows selecting genes for diff analysis based upon minimum counts across
#' samples within a group (e.g. minimum counts across replicate samples).
#'
#' @param counts a character string providing a name of count file or a data frame of counts for each gene.
#' See details below to know more about format of file and count data frame.
#' @param column_geneid a character string denoting a column of geneid in \code{counts}
#' @param column_samples a character vector denoting names of sample columns from \code{counts}
#' @param sample_info a character string denoting a name of sample information file or a data frame.
#' A file or a data frame both must have at least two columns without column names. First column denotes to samples names
#' and second column denotes group name for each sample in first column. For e.g.
#'
#'  |       |  |
#'  | ----------- | ----------- |
#'  | sample 1      | WT       |
#'  | sample 2      | WT       |
#'  | sample 3      | KO       |
#'  | sample 4      | KO       |
#'
#' @param group_numerator a character vector denoting sample groups to use in numerator to calculate fold change.
#' @param group_denominator a character vector denoting sample groups to use in denominator to calculate fold change.
#' @param delim a character denoting deliminator for `count` file. Only valid if `count` is a file path.
#' @param comment_char a character denoting comments line in count file. Only valid if `count` is a file path.
#' @param min_counts a numeric value, default 10,  denoting minimum counts for a gene to be used o consider a
#' gene for differential expression analysis.
#' @param min_replicates a numeric value, default 1, denoting minimum samples within a group must have `minimum_counts`.
#' Value provided must not be higher than number of replicates in each group.
#' For example for given values `min_replicates = 2` and  `minimum_counts = 10`
#' the genes which have minimum counts 10 in atleast 2 sample groups will be used for DEG.
#' @param cutoff_lfc a numeric value which is internally passed to \link{categorize_diff_genes}
#' @param cutoff_padj a numeric value which is internally passed to \link{categorize_diff_genes}
#' @param ... for future use
#' @param cutoff_pval a numeric value which is internally passed to \link{categorize_diff_genes}
#' @param regul_based_upon either of 1, 2 or 3 which is internally passed to \link{categorize_diff_genes}
#'
#' @return a data frame of DESeq results, DEG, and DEG summary.
#' @export
#' @importFrom cli cli_alert
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_alert_success
#' @importFrom cli cli_alert_warning
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 results
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate_if
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr tally
#' @importFrom dplyr ungroup
#' @importFrom glue glue
#' @importFrom magrittr `%<>%`
#' @importFrom magrittr `%>%`
#' @importFrom purrr map
#' @importFrom readr read_delim
#' @importFrom rlang `:=`
#' @importFrom rlang `!!!`
#' @importFrom rlang `!!`
#' @importFrom rlang enquo
#' @importFrom rlang enquos
#' @importFrom rlang is_scalar_double
#' @importFrom rlang new_formula
#' @importFrom rlang quo
#' @importFrom rlang sym
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr expand_grid
#' @importFrom tidyr gather
#' @import DESeq2
#' @import TidyWrappers
#' @examples
#' \dontrun{
#'
#' set.seed(123)
#' # create dummy RNAseq (dummy) count matrix
#' counts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) %>%
#'   as.data.frame() %>% tibble::as_tibble()
#'
#' colnames(counts) <- c(paste("c" , c(1:5), sep = ""),c(paste ("d" , 1:5, sep = "")))
#'
#' counts %<>%  dplyr::mutate("Geneid" = stringi::stri_rand_strings(n = 100, length = 5))  %>%
#'   dplyr::relocate("Geneid")
#' # create sample info
#' si <- tibble::tibble(samples = colnames(counts)[-1] , sample_groups = factor(rep(c("c","d"), each=5)))
#'
#' res <- get_deg(counts = counts ,
#'                sample_info = si,
#'                column_geneid = "Geneid" , group_numerator = "d" , group_denominator = "c",
#'                column_samples = c("c1","c2","c3","c4" ,"c5" ,"d1" ,"d2","d3" ,"d4" ,"d5"))
#'
#' names <- paste(res$)
#'
#' ## DESEq result object(s)
#' res$dsr
#'
#' ## DESEq result data frame
#' res$dsr_tibble
#'
#' ## DESEq result data frame  DEG assigned, look at the columns 'signif' and 'regul'
#'
#' res$dsr_tibble_deg
#'
#' res$dsr_tibble_deg[[1]] %>% dplyr::filter(regul != "other")
#'
#' ## DEG summary
#'
#' res$deg_summmary
#'
#'
#' }
get_deg <- function(counts,
                    column_geneid,
                    column_samples,
                    sample_info,
                    group_numerator,
                    group_denominator,
                    delim = "\t", # delim for count file. valid only if count file provided. default "\t".
                    comment_char = "#", # comment char for count file. valid only if count file provided. default "#",
                    min_counts = 10, # minimum counts in anyone sample
                    min_replicates = 1,
                    cutoff_lfc = 1,
                    cutoff_pval = 0.05,
                    cutoff_padj = 0.01,
                    regul_based_upon = 1,
                    ...
                    ){

  stop("`get_deg` is deprecated.
             Please use `run_deseq_analysis` instead.")


  # defaults
  # counts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) %>% as.data.frame() %>% tibble::as_tibble()

  # colnames(counts) <- c(paste("c" , c(1:5), sep = ""),c(paste ("d" , 1:5, sep = "")))

  # counts %<>%  dplyr::mutate("Geneid" = stringi::stri_rand_strings(n = 100, length = 5))  %>% dplyr::relocate("Geneid")
  # sample_info <- tibble::tibble(samples = colnames(counts)[-1] , sample_groups = factor(rep(c("c","d"), each=5)))
  # column_geneid = "Geneid"
  # column_samples = c("c1","c2","c3","c4" ,"c5" ,"d1" ,"d2","d3" ,"d4" ,"d5")
  # cutoff_pval = 0.05
  # cutoff_padj = 0.01
  # cutoff_lfc = 1.5
  # delim = "\t"
  # comment_char = "#"
  # min_counts <- 10 # minimum counts in anyone sample
  # min_replicates <- 2 # number of replicates fulfilling criteria of minimum counts in anyone sample (min_counts)
  # group_numerator = "d"
  # group_denominator = "c"

  ## define statics/global
  sample_info_colnames = c("sample_names", "sample_groups")

  # validate numeric and scaler value
  num_args = list(min_counts, min_replicates, cutoff_lfc , cutoff_padj)
  num_args_quo <- rlang::quos(min_counts, min_replicates, cutoff_lfc , cutoff_padj)
  names(num_args)  <- purrr::map_chr(num_args_quo , ~ rlang::quo_name(..1))

  is_numeric_and_scaler <- purrr::imap_lgl(num_args ,~dplyr::if_else(rlang::is_scalar_double(..1), TRUE, FALSE))

  for(i in seq_along(is_numeric_and_scaler)){
    if(!is_numeric_and_scaler[i])
      stop(glue::glue("{names(is_numeric_and_scaler[i])} must be a double of length 1."))
  }

  if(!rlang::is_scalar_double(min_counts)){
    stop("'min_counts' must be of length one of class double.")
  }


  ## prepare count matrix
  if(!inherits(counts , what = c("character","tbl_df","tbl" ,"data.frame"))){
    stop("counts must be an object of one of these classes - character, tbl_df, tbl, or data.frame")
  }

  if(inherits(counts , what = c("character"))){
    cli::cli_alert_info(" Reading count file ...")
    count_data <- readr::read_delim(counts ,
                                    delim = delim,
                                    comment = comment_char,
                                    trim_ws = TRUE,
                                    col_names = TRUE,
                                    show_col_types = F ,
                                    progress = T)
    cli::cli_alert_success("Done.")
  } else {
    count_data <- counts
  }

  # remove duplicate rows if any
  count_data %<>% unique()

  # value for an argument column_geneid must be a character string

  if(!is.character(column_geneid) && length(column_geneid) != 1){
    stop("Value for an argument 'geneid_coulmn' must be a character string")
  }

  # value for an argument columns_samples must be a character vector

  if(!is.character(column_samples) && length(column_samples) > 1){
    stop("Value for an argument 'sample_columns' must be a character vector of length > 1")
  }

  # check if columns specified in column_samples are present in count_data
  count_data_col_names <-  colnames(count_data)

  if(!all(column_samples %in% count_data_col_names)){
    index_not_present <- which(!(column_samples %in% count_data_col_names))
    value_not_present <- column_samples[index_not_present] %>% paste0(collapse = "','")
    stop(glue::glue("Values '{value_not_present}' in 'column_samples' must present as column(s) in 'counts'."))
  }

  # quote column names
  column_geneid_quo <- rlang::enquo(column_geneid)
  column_samples_quo <- rlang::enquos(column_samples)

  ## prepare sample information

  # value for sample_info either from a class character, tbl_df, tbl or data.frame

  if(!inherits(sample_info, what = c("character","tbl_df","tbl" ,"data.frame"))){
    stop("sample_info must be an object of one of these classes - character, tbl_df, tbl, or data.frame")
  }

  # if sample_info is character then read data from file.

  if(inherits(sample_info , what = c("character"))){
    cli::cli_alert_info(" Reading sample_info file ...")
    sample_info_data <- readr::read_delim(sample_info ,
                                          delim = "\t",
                                          trim_ws = TRUE,
                                          col_names = FALSE,
                                          col_types = c("cc"), # both column will be considered as character
                                          show_col_types = FALSE,
                                          progress = TRUE)

    cli::cli_alert_success("Done.")
  } else{
    sample_info_data <- sample_info
  }
  # if more than 2 columns present in sample_info keep only 2 and show warning.
  if(ncol(sample_info_data) != 2 ){
    if(ncol(sample_info_data) > 2 ){
      cli::cli_alert_warning("'sample_info' contains more than 2 columns. Only first 2 will be considered.")
      sample_info_data %<>%
        dplyr::select(1:2)
    } else if(ncol(sample_info_data) < 2 ) {
      stop("'sample_info' must contains 2 columns deliminated by '\t'")
    }
  } else{
    sample_info_data = sample_info_data
  }

  # remove duplicates and keep only first two columns.
  sample_info_data %<>%
    dplyr::select(1:2) %>% # if more than 2 columns, keep only first 2.
    unique() # remove duplicate rows if any

  # assign column names to sample_info_data
  colnames(sample_info_data) <- sample_info_colnames

  # check if all the samples provided in sample_info are present in count_data file

  sample_info_samples <- sample_info_data %>% dplyr::pull(sample_info_colnames[1])
  sample_info_samples_quo <- rlang::enquos(sample_info_samples)

  if(!all(sample_info_samples %in% count_data_col_names)){
    index_not_present <- which(!(sample_info_samples %in% count_data_col_names))
    value_not_present <- sample_info_samples[index_not_present] %>% paste0(collapse = "','")
    stop(glue::glue("Values '{value_not_present}' given in first column of 'sample_info' must present as column(s) in 'counts'."))
  }

  # validate 'group_numerator' and 'group_denominator'

  if(!inherits(group_numerator , what = "character")) {
    stop("'group_numerator' must be a character vector.")
  }

  if(!inherits(group_denominator , what = "character")) {
    stop("'group_denominator' must be a character vector.")
  }

  # values mentioned for 'group_numerator' and 'group_denominator' must present in sample_info
  sample_info_groups <- sample_info_data %>% dplyr::pull(sample_info_colnames[2]) %>% unique()

  if(!all(group_numerator %in% sample_info_groups)){
    stop("all values from 'group_numerator' must present in 'sample_info' column 2.")
  }

  if(!all(group_denominator %in% sample_info_groups)){
    stop("all values from 'group_denominator' must present in 'sample_info' column 2.")
  }


  # subset necessary columns from count.
  #The samples which are present in sample_information will only be processed further.

  count_data_select <- count_data %>%
    dplyr::select(!!column_geneid_quo, !!!sample_info_samples_quo)

  # display warning for the samples which won't be considered for analysis

  # // TO DO

  # in count_data convert NA to 0

  count_data_select_NA <- count_data_select %>% TidyWrappers::tbl_keep_rows_NA_any()
  if(nrow(count_data_select_NA) >= 1) {
    cli::cli_alert_warning("Following gene(s) contains NA in one or more samples. Either remove tem or Default NA will be converted to 0.")
    genes_NA <- count_data_select_NA %>%
      dplyr::pull(!!column_geneid_quo)
    cli::cli_alert(genes_NA)
    # replace NA with 0
    count_data_select %<>% dplyr::mutate_if(is.numeric , ~ ..1 %>% tidyr::replace_na(0))
  }

  ## prepare annotation
  ## // TO DO
  # cli::cli_alert_info(text = "Preparing annotations")
  # features_slected <- filter_gff(gtf_file = gtf_file)
  # cli::cli_alert_success(text = "Done")

  # Run DESeq on those genes which fits on user defined cutoffs - minimum counts to each gene and minimum number of replicates with minimum counts

  count_data_all_zero <- count_data_select %>%
    TidyWrappers::tbl_keep_rows_zero_all()

  count_data_non_zero_in_one <- count_data_select %>%
    TidyWrappers::tbl_remove_rows_zero_all()

  # to perform DESeq analysis genes having non zero count in at least one sample will be used.

  genes_to_keep <- count_data_non_zero_in_one %>%

    # convert data wide to long format.
    tidyr::gather(-!!column_geneid_quo, key = "sample", value = "counts") %>%

    # join sampleinfo table
    dplyr::left_join(sample_info_data, by = c("sample" = sample_info_colnames[1])) %>%

    # to identify genes which has counts more than cutoff across replicates group by sample groups followed by gene id.
    # total groups will be number of unique genes * number of sample groups.
    dplyr::group_by(!!rlang::sym(sample_info_colnames[2]), !!rlang::sym(column_geneid)) %>%

    dplyr::arrange(!!rlang::sym(column_geneid), !!rlang::sym(sample_info_colnames[2])) %>%

    # add column (count_above_threshold) containing counts of replicates having min counts >= min_counts
    dplyr::mutate(count_above_threshold = sum(counts >= min_counts)) %>%

    # filter genes count_above_thresholds >= min_replicates (e.g. 2);
    # meaning that for a given gene at least two replicates have count >= threshold (e.g., 10)
    dplyr::filter(count_above_threshold >= min_replicates) %>%
    dplyr::ungroup() %>%
    dplyr::pull(!!rlang::sym(column_geneid)) %>%
    unique()

  #genes_to_keep %>% length()

  count_data_filter <- count_data_select %>%
    dplyr::filter(!!rlang::sym(column_geneid) %in% genes_to_keep)

  ## DESeq 2 analysis
  # in order to run DESeq2 order of rows in sample_info_data and columns in count must be same.

  if(!all(sample_info_data[[1]] == colnames(count_data_filter)[-1])){
    index <- match(sample_info_data[[1]], colnames(count_data_filter))
    count_data_filter %<>% dplyr::select(!!column_geneid,dplyr::all_of(index))
  }

  # run DEseq
  cli::cli_alert_info("Running DESeq2 ...")

  # create DESeqDataSet object -
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data_filter %>% as.data.frame(),
                                        colData = sample_info_data %>% tibble::column_to_rownames(sample_info_colnames[[1]]),
                                        design =  rlang::new_formula(NULL, rlang::sym(sample_info_colnames[2])) ,
                                        tidy = T)

  dds <- DESeq2::DESeq(dds)

  cli::cli_alert_success("Done.")

  # generate normalize count matrix

  norm_counts <- DESeq2::counts(dds, normalized  = T) %>%  as.data.frame() %>%
    tibble::rownames_to_column(var = column_geneid) %>% tibble::as_tibble()

  # create combinations of all comparisons from  group_numerator and group_denominator.


  comb <- tidyr::expand_grid(x = group_numerator, y = group_denominator)

  column_numerator = rlang::quo(numerator)
  column_denominator = rlang::quo(denominator)
  comb %<>% dplyr::rename(!!column_numerator := x , !!column_denominator := y )

  # remove all combinations where numerator and denominator
  comb <- comb %>% dplyr::filter(numerator != denominator)

  if(nrow(comb) == 0){
    stop("Numerator and Denominator have at least a pair of different value.")
  }

  # summarize results

  cli::cli_alert_info("Summarizing DEG ...")

  # prepare a tibble where each row denotes to a combination of numerator and denominator.

  xx <- comb %>%

    # for each combination of numerator and denominator get deseq result. Results will be stored in a list column of tibble

    dplyr::mutate(dsr =  purrr::map2( .x = !!column_numerator,
                                      .y = !!column_denominator ,
                                      ~ DESeq2::results(dds, contrast = c("sample_groups" , ..1,..2)))) %>%

    # convert deseq result object into a tibble.

    dplyr::mutate(dsr_tibble = purrr::map(dsr , ~..1 %>% dsr_to_tibble())) %>%

    # categorize DEG based on log2FC and pvalue

    dplyr::mutate(dsr_tibble_deg = purrr::map(dsr_tibble , ~ ..1 %>%
                                                categorize_diff_genes(log2fc_cutoff = cutoff_lfc ,
                                                                      pval_cutoff = cutoff_pval,
                                                                      padj_cutoff = cutoff_padj ,
                                                                      regul_based_upon = regul_based_upon) )) %>%
    # add count matrix

    dplyr::mutate(norm_counts = list(norm_counts = norm_counts)) %>%

    # summarize DEG

    dplyr::mutate(deg_summmary = purrr::map(dsr_tibble_deg , ~ ..1 %>%
                                              dplyr::group_by(regul) %>%
                                              dplyr::tally() ,.id = "cond"))



  cli::cli_alert_info("Done.")
    return(xx)

}



#' @title Perform differential expression analysis using [DESeq2::DESeq()]
#' @description This is a wrapper function build upon [DESeq2::DESeq()]  and  [DESeq2::DESeqResults()]
#' to find diff genes and categories them based on various cutoffs such p-value, padj-value, log2fc etc.
#' It also allows selecting genes for diff analysis based upon minimum counts across
#' samples within a group (e.g. minimum counts across replicate samples).
#'
#' @param counts a character string providing a name of count file or a data frame of counts for each gene.
#' See details below to know more about format of file and count data frame.
#' @param column_geneid a character string denoting a column of geneid in \code{counts}
#' @param column_samples a character vector denoting names of sample columns from \code{counts}
#' @param sample_info a character string denoting a name of sample information file or a data frame.
#' A file or a data frame both must have at least two columns without column names. First column denotes to samples names
#' and second column denotes group name for each sample in first column. For e.g.
#'
#'  |       |  |
#'  | ----------- | ----------- |
#'  | sample 1      | WT       |
#'  | sample 2      | WT       |
#'  | sample 3      | KO       |
#'  | sample 4      | KO       |
#'
#' @param group_numerator a character vector denoting sample groups to use in numerator to calculate fold change.
#' @param group_denominator a character vector denoting sample groups to use in denominator to calculate fold change.
#' @param delim a character denoting deliminator for `count` file. Only valid if `count` is a file path.
#' @param comment_char a character denoting comments line in count file. Only valid if `count` is a file path.
#' @param min_counts a numeric value, default 10,  denoting minimum counts for a gene to be used o consider a
#' gene for differential expression analysis.
#' @param min_replicates a numeric value, default 1, denoting minimum samples within a group must have `minimum_counts`.
#' Value provided must not be higher than number of replicates in each group.
#' For example for given values `min_replicates = 2` and  `minimum_counts = 10`
#' the genes which have minimum counts 10 in atleast 2 sample groups will be used for DEG.
#' @param cutoff_lfc a numeric value which is internally passed to \link{categorize_diff_genes}
#' @param cutoff_padj a numeric value which is internally passed to \link{categorize_diff_genes}
#' @param ... for future use
#' @param cutoff_pval a numeric value which is internally passed to \link{categorize_diff_genes}
#' @param regul_based_upon either of 1, 2 or 3 which is internally passed to \link{categorize_diff_genes}
#'
#' @return a data frame of DESeq results, DEG, and DEG summary.
#' @export
#' @importFrom cli cli_alert
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_alert_success
#' @importFrom cli cli_alert_warning
#' @importFrom DESeq2 DESeq
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 results
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate_if
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr tally
#' @importFrom dplyr ungroup
#' @importFrom glue glue
#' @importFrom magrittr `%<>%`
#' @importFrom magrittr `%>%`
#' @importFrom purrr map
#' @importFrom readr read_delim
#' @importFrom rlang `:=`
#' @importFrom rlang `!!!`
#' @importFrom rlang `!!`
#' @importFrom rlang enquo
#' @importFrom rlang enquos
#' @importFrom rlang is_scalar_double
#' @importFrom rlang new_formula
#' @importFrom rlang quo
#' @importFrom rlang sym
#' @importFrom tibble column_to_rownames
#' @importFrom tidyr expand_grid
#' @importFrom tidyr gather
#' @import DESeq2
#' @import TidyWrappers
#' @examples
#' \dontrun{
#'
#' set.seed(123)
#' # create dummy RNAseq (dummy) count matrix
#' counts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) %>%
#'   as.data.frame() %>% tibble::as_tibble()
#'
#' colnames(counts) <- c(paste("c" , c(1:5), sep = ""),c(paste ("d" , 1:5, sep = "")))
#'
#' counts %<>%  dplyr::mutate("Geneid" = stringi::stri_rand_strings(n = 100, length = 5))  %>%
#'   dplyr::relocate("Geneid")
#' # create sample info
#' si <- tibble::tibble(samples = colnames(counts)[-1] , sample_groups = factor(rep(c("c","d"), each=5)))
#'
#' res <- run_deseq_analysis(counts = counts ,
#'                sample_info = si,
#'                column_geneid = "Geneid" , group_numerator = "d" , group_denominator = "c",
#'                column_samples = c("c1","c2","c3","c4" ,"c5" ,"d1" ,"d2","d3" ,"d4" ,"d5"))
#'
#' names <- paste(res$)
#'
#' ## DESEq result object(s)
#' res$dsr
#'
#' ## DESEq result data frame
#' res$dsr_tibble
#'
#' ## DESEq result data frame  DEG assigned, look at the columns 'signif' and 'regul'
#'
#' res$dsr_tibble_deg
#'
#' res$dsr_tibble_deg[[1]] %>% dplyr::filter(regul != "other")
#'
#' ## DEG summary
#'
#' res$deg_summmary
#'
#'}
#'
run_deseq_analysis <- function(
  counts,
  column_geneid,
  column_samples,
  sample_info,
  group_numerator,
  group_denominator,
  delim = "\t", # delim for count file. valid only if count file provided. default "\t".
  comment_char = "#", # comment char for count file. valid only if count file provided. default "#",
  min_counts = 10, # minimum counts in anyone sample
  min_replicates = 1,
  cutoff_lfc = 1,
  cutoff_pval = 0.05,
  cutoff_padj = 0.01,
  regul_based_upon = 1,
  ...

){

  # defaults
   # counts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) %>% as.data.frame() %>% tibble::as_tibble()
   #
   # colnames(counts) <- c(paste("c" , c(1:5), sep = ""),c(paste ("d" , 1:5, sep = "")))
   #
   # counts %<>%  dplyr::mutate("Geneid" = stringi::stri_rand_strings(n = 100, length = 5))  %>% dplyr::relocate("Geneid")
   # sample_info <- tibble::tibble(samples = colnames(counts)[-1] , sample_groups = factor(rep(c("c","d"), each=5)))
   # column_geneid = "Geneid"
   # column_samples = c("c1","c2","c3","c4" ,"c5" ,"d1" ,"d2","d3" ,"d4" ,"d5")
   # cutoff_pval = 0.05
   # cutoff_padj = 0.01
   # cutoff_lfc = 1.5
   # delim = "\t"
   # comment_char = "#"
   # min_counts <- 10 # minimum counts in anyone sample
   # min_replicates <- 2 # number of replicates fulfilling criteria of minimum counts in anyone sample (min_counts)
   # group_numerator = "d"
   # group_denominator = "c"

  ## define statics/global
  sample_info_colnames = c("sample_names", "sample_groups")

  # validate numeric and scaler value
  num_args = list(min_counts, min_replicates, cutoff_lfc , cutoff_padj)
  num_args_quo <- rlang::quos(min_counts, min_replicates, cutoff_lfc , cutoff_padj)
  names(num_args)  <- purrr::map_chr(num_args_quo , ~ rlang::quo_name(..1))

  is_numeric_and_scaler <- purrr::imap_lgl(num_args ,~dplyr::if_else(rlang::is_scalar_double(..1), TRUE, FALSE))

  for(i in seq_along(is_numeric_and_scaler)){
    if(!is_numeric_and_scaler[i])
      stop(glue::glue("{names(is_numeric_and_scaler[i])} must be a double of length 1."))
  }

  if(!rlang::is_scalar_double(min_counts)){
    stop("'min_counts' must be of length one of class double.")
  }


  ## prepare count matrix
  if(!inherits(counts , what = c("character","tbl_df","tbl" ,"data.frame"))){
    stop("counts must be an object of one of these classes - character, tbl_df, tbl, or data.frame")
  }

  if(inherits(counts , what = c("character"))){
    cli::cli_alert_info(" Reading count file ...")
    count_data <- readr::read_delim(counts ,
                                    delim = delim,
                                    comment = comment_char,
                                    trim_ws = TRUE,
                                    col_names = TRUE,
                                    show_col_types = F ,
                                    progress = T)
    cli::cli_alert_success("Done.")
  } else {
    count_data <- counts
  }

  # remove duplicate rows if any
  count_data %<>% unique()

  # value for an argument column_geneid must be a character string

  if(!is.character(column_geneid) && length(column_geneid) != 1){
    stop("Value for an argument 'geneid_coulmn' must be a character string")
  }

  # value for an argument columns_samples must be a character vector

  if(!is.character(column_samples) && length(column_samples) > 1){
    stop("Value for an argument 'sample_columns' must be a character vector of length > 1")
  }

  # check if columns specified in column_samples are present in count_data
  count_data_col_names <-  colnames(count_data)

  if(!all(column_samples %in% count_data_col_names)){
    index_not_present <- which(!(column_samples %in% count_data_col_names))
    value_not_present <- column_samples[index_not_present] %>% paste0(collapse = "','")
    stop(glue::glue("Values '{value_not_present}' in 'column_samples' must present as column(s) in 'counts'."))
  }

  # quote column names
  column_geneid_quo <- rlang::enquo(column_geneid)
  column_samples_quo <- rlang::enquos(column_samples)

  ## prepare sample information

  # value for sample_info either from a class character, tbl_df, tbl or data.frame

  if(!inherits(sample_info, what = c("character","tbl_df","tbl" ,"data.frame"))){
    stop("sample_info must be an object of one of these classes - character, tbl_df, tbl, or data.frame")
  }

  # if sample_info is character then read data from file.

  if(inherits(sample_info , what = c("character"))){
    cli::cli_alert_info(" Reading sample_info file ...")
    sample_info_data <- readr::read_delim(sample_info ,
                                          delim = "\t",
                                          trim_ws = TRUE,
                                          col_names = FALSE,
                                          col_types = c("cc"), # both column will be considered as character
                                          show_col_types = FALSE,
                                          progress = TRUE)

    cli::cli_alert_success("Done.")
  } else{
    sample_info_data <- sample_info
  }
  # if more than 2 columns present in sample_info keep only 2 and show warning.
  if(ncol(sample_info_data) != 2 ){
    if(ncol(sample_info_data) > 2 ){
      cli::cli_alert_warning("'sample_info' contains more than 2 columns. Only first 2 will be considered.")
      sample_info_data %<>%
        dplyr::select(1:2)
    } else if(ncol(sample_info_data) < 2 ) {
      stop("'sample_info' must contains 2 columns deliminated by '\t'")
    }
  } else{
    sample_info_data = sample_info_data
  }

  # remove duplicates and keep only first two columns.
  sample_info_data %<>%
    dplyr::select(1:2) %>% # if more than 2 columns, keep only first 2.
    unique() # remove duplicate rows if any

  # assign column names to sample_info_data
  colnames(sample_info_data) <- sample_info_colnames

  # check if all the samples provided in sample_info are present in count_data file

  sample_info_samples <- sample_info_data %>% dplyr::pull(sample_info_colnames[1])
  sample_info_samples_quo <- rlang::enquos(sample_info_samples)

  if(!all(sample_info_samples %in% count_data_col_names)){
    index_not_present <- which(!(sample_info_samples %in% count_data_col_names))
    value_not_present <- sample_info_samples[index_not_present] %>% paste0(collapse = "','")
    stop(glue::glue("Values '{value_not_present}' given in first column of 'sample_info' must present as column(s) in 'counts'."))
  }

  # validate 'group_numerator' and 'group_denominator'

  if(!inherits(group_numerator , what = "character")) {
    stop("'group_numerator' must be a character vector.")
  }

  if(!inherits(group_denominator , what = "character")) {
    stop("'group_denominator' must be a character vector.")
  }

  # values mentioned for 'group_numerator' and 'group_denominator' must present in sample_info
  sample_info_groups <- sample_info_data %>% dplyr::pull(sample_info_colnames[2]) %>% unique()

  if(!all(group_numerator %in% sample_info_groups)){
    stop("all values from 'group_numerator' must present in 'sample_info' column 2.")
  }

  if(!all(group_denominator %in% sample_info_groups)){
    stop("all values from 'group_denominator' must present in 'sample_info' column 2.")
  }


  # subset necessary columns from count.
  #The samples which are present in sample_information will only be processed further.

  count_data_select <- count_data %>%
    dplyr::select(!!column_geneid_quo, !!!sample_info_samples_quo)

  # display warning for the samples which won't be considered for analysis

  # // TO DO

  # in count_data convert NA to 0

  count_data_select_NA <- count_data_select %>% TidyWrappers::tbl_keep_rows_NA_any()
  if(nrow(count_data_select_NA) >= 1) {
    cli::cli_alert_warning("Following gene(s) contains NA in one or more samples. Either remove tem or Default NA will be converted to 0.")
    genes_NA <- count_data_select_NA %>%
      dplyr::pull(!!column_geneid_quo)
    cli::cli_alert(genes_NA)
    # replace NA with 0
    count_data_select %<>% dplyr::mutate_if(is.numeric , ~ ..1 %>% tidyr::replace_na(0))
  }

  ## prepare annotation
  ## // TO DO
  # cli::cli_alert_info(text = "Preparing annotations")
  # features_slected <- filter_gff(gtf_file = gtf_file)
  # cli::cli_alert_success(text = "Done")

  # Run DESeq on those genes which fits on user defined cutoffs - minimum counts to each gene and minimum number of replicates with minimum counts

  count_data_all_zero <- count_data_select %>%
    TidyWrappers::tbl_keep_rows_zero_all()

  count_data_non_zero_in_one <- count_data_select %>%
    TidyWrappers::tbl_remove_rows_zero_all()

  # to perform DESeq analysis genes having non zero count in at least one sample will be used.

  genes_to_keep <- count_data_non_zero_in_one %>%

    # convert data wide to long format.
    tidyr::gather(-!!column_geneid_quo, key = "sample", value = "counts") %>%

    # join sampleinfo table
    dplyr::left_join(sample_info_data, by = c("sample" = sample_info_colnames[1])) %>%

    # to identify genes which has counts more than cutoff across replicates group by sample groups followed by gene id.
    # total groups will be number of unique genes * number of sample groups.
    dplyr::group_by(!!rlang::sym(sample_info_colnames[2]), !!rlang::sym(column_geneid)) %>%

    dplyr::arrange(!!rlang::sym(column_geneid), !!rlang::sym(sample_info_colnames[2])) %>%

    # add column (count_above_threshold) containing counts of replicates having min counts >= min_counts
    dplyr::mutate(count_above_threshold = sum(counts >= min_counts)) %>%

    # filter genes count_above_thresholds >= min_replicates (e.g. 2);
    # meaning that for a given gene at least two replicates have count >= threshold (e.g., 10)
    dplyr::filter(count_above_threshold >= min_replicates) %>%
    dplyr::ungroup() %>%
    dplyr::pull(!!rlang::sym(column_geneid)) %>%
    unique()

  #genes_to_keep %>% length()

  count_data_filter <- count_data_select %>%
    dplyr::filter(!!rlang::sym(column_geneid) %in% genes_to_keep)

  ## DESeq 2 analysis
  # in order to run DESeq2 order of rows in sample_info_data and columns in count must be same.

  if(!all(sample_info_data[[1]] == colnames(count_data_filter)[-1])){
    index <- match(sample_info_data[[1]], colnames(count_data_filter))
    count_data_filter %<>% dplyr::select(!!column_geneid,dplyr::all_of(index))
  }

  # run DEseq
  cli::cli_alert_info("Running DESeq2 ...")

  # create DESeqDataSet object -
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_data_filter %>% as.data.frame(),
                                        colData = sample_info_data %>% tibble::column_to_rownames(sample_info_colnames[[1]]),
                                        design =  rlang::new_formula(NULL, rlang::sym(sample_info_colnames[2])) ,
                                        tidy = T)

  dds <- DESeq2::DESeq(dds)

  cli::cli_alert_success("Done.")

  # generate normalize count matrix

  norm_counts <- DESeq2::counts(dds, normalized  = T) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = column_geneid) %>%
    tibble::as_tibble()

  # group normalized counts for each sample

  norm_counts <- norm_counts %>%
    tidyr::pivot_longer(cols = -!!column_geneid , names_to = "samples", values_to = "norm_counts") %>%
    dplyr::left_join(sample_info_data , by = c("samples" = "sample_names"))  %>%
    parcutils::named_group_split(sample_groups) %>%
    purrr::map(~ ..1 %>% tidyr::pivot_wider(id_cols = !!column_geneid,
                                            names_from = samples ,
                                            values_from = norm_counts))

  # create combinations of all comparisons from  group_numerator and group_denominator.


  comb <- tidyr::expand_grid(x = group_numerator, y = group_denominator)

  column_numerator = rlang::quo(numerator)
  column_denominator = rlang::quo(denominator)
  comb %<>% dplyr::rename(!!column_numerator := x , !!column_denominator := y )

  # remove all combinations where numerator and denominator
  comb <- comb %>% dplyr::filter(numerator != denominator)

  if(nrow(comb) == 0){
    stop("Numerator and Denominator have at least a pair of different value.")
  }

  # summarize results

  cli::cli_alert_info("Summarizing DEG ...")

  # prepare a tibble where each row denotes to a combination of numerator and denominator.

  xx <- comb %>%

    ## add column comparisons
    dplyr::mutate(comp = stringr::str_c(.$numerator , .$denominator ,sep = "_VS_")) %>%

    ## add column norm counts
    dplyr::mutate(norm_counts = purrr::map2(.$numerator, .$denominator, ~ norm_counts[c(..1,..2)])) %>%


    # for each combination of numerator and denominator get deseq result. Results will be stored in a list column of tibble

    dplyr::mutate(dsr =  purrr::map2( .x = !!column_numerator,
                                      .y = !!column_denominator ,
                                      ~ DESeq2::results(dds, contrast = c("sample_groups" , ..1,..2)))) %>%

    # convert deseq result object into a tibble.

    dplyr::mutate(dsr_tibble = purrr::map(dsr , ~..1 %>% dsr_to_tibble())) %>%

    # categorize DEG based on log2FC and pvalue

    dplyr::mutate(dsr_tibble_deg = purrr::map(dsr_tibble , ~ ..1 %>%
                                                categorize_diff_genes(log2fc_cutoff = cutoff_lfc ,
                                                                      pval_cutoff = cutoff_pval,
                                                                      padj_cutoff = cutoff_padj ,
                                                                      regul_based_upon = regul_based_upon) )) %>%
    # summarize DEG

    dplyr::mutate(deg_summmary = purrr::map(dsr_tibble_deg , ~ ..1 %>%
                                              dplyr::group_by(regul) %>%
                                              dplyr::tally() ,.id = "cond"))


  ## use comparisons as names for list elements
  xx  %<>%  dplyr::mutate_if(is.list, ~rlang::set_names(. , xx$comp))

  # assign names to each list in the tibble

  cli::cli_alert_info("Done.")

  ## set class

  class(xx) <- c("parcutils" , class(xx))
  return(xx)

}



#' Convert DESeq result object in to a tibble.
#'
#' @param x [DESeq2::DESeqResults()] object
#'
#' @return a data frame
#' @export
#' @keywords internal
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate_if
#' @importFrom tibble as_tibble
#' @examples
#' \dontrun{
#'  // TO DO
#' }
#'
dsr_to_tibble <- function(x, .col_gene_id = "gene_id"){

  x %>% as.data.frame() %>%
    tibble::rownames_to_column(.col_gene_id) %>%
    tibble::as_tibble() #%>%
    #dplyr::mutate_if(is.numeric, ~ round(..1, 3)) ## round by 3 digits
}



##@@@@@@@@
## filter gff file
##@@@@@@@@

#' Filter rows from GFF file.
#'
#' @param gtf_file  // TO DO
#' @param feature_type // TO DO
#' @param feature_biotype // TO DO
#' @param column_gene_id // TO DO
#' @param column_feature_type // TO DO
#' @param column_feature_biotype // TO DO
#'
#' @return // TO DO
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' // TO DO
#' }
#'
filter_gff <- function(gtf_file ,
                       feature_type = "gene" ,
                       feature_biotype = "protein_coding",
                       column_gene_id = "gene_id",
                       column_feature_type = "type" ,
                       column_feature_biotype = "gene_biotype"){

  # gtf_file = gtf_file
  # feature_type = feature_type
  # feature_bio_type = feature_bio_type

  column_feature_type_quo = rlang::enquo(column_feature_type)
  column_feature_biotype_quo = rlang::enquo(column_feature_biotype)
  column_gene_id_quo  = rlang::enquo(column_gene_id)

  cli::cli_alert_info(text = "Reading gft/gff file ...")
  gtf_annot <- rtracklayer::import(gtf_file) %>%
    as.data.frame() %>%
    tibble::as_tibble()
  cli::cli_alert_success(text = "Done.")

  # get genes
  gtf_annot_subset <- gtf_annot %>%
    dplyr::filter(type == feature_type) %>%
    dplyr::filter(gene_biotype == feature_bio_type) %>%
    dplyr::select(!!column_gene_id_quo, !!column_feature_type_quo, !!column_feature_biotype_quo)

  return(gtf_annot_subset)
}


#' @title  Categorize diff genes based on log2FC, p-value and p-adj value.
#' @description // TO DO
#'
#' @param dsr_tibble  a data frame obtained from DESeqResult object or an object of DESeqResult.
#' @param log2fc_cutoff,pval_cutoff,padj_cutoff a numeric value, default `(log2fc_cutoff = 1, pval_cutoff =  0.05, padj_cutoff = 0.01)`
#' denoting cutoffs. These criteria will be used to decide significance  `(NS, p-value, log2FC, p-value&log2FC)` and
#' type of regulation `(Up, Down & other)` of diff genes.
#' See details.
#' @param add_column_regul logical, default `TRUE`, indicating whether to add column  `add_column_regul` or not. Added column will contain the values
#' Up', 'Down' or 'other'.
#' @param regul_based_upon one of the numeric choices  1, 2, or 3.
#' ## if 1 ...
#'  - Up : log2fc >= log2fc_cutoff & pvalue <= pvalue_cutoff
#'  - Down : log2fc  <= (-1) * log2fc_cutoff & pvalue <= pvalue_cutoff
#'  - Other : remaining genes
#' ## if 2 ...
#'  - Up : log2fc >= log2fc_cutoff & padj <= padj_cutoff
#'  - Down : log2fc  <= (-1) * log2fc_cutoff & padj <= padj_cutoff
#'  - Other : remaining genes
#' ## if 3 ...
#'  - Up : log2fc >= log2fc_cutoff & pvalue <= pvalue_cutoff & padj <= padj_cutoff
#'  - Down : log2fc  <= (-1) * log2fc_cutoff pvalue <= pvalue_cutoff & padj <= padj_cutoff
#'  - Other : remaining genes
#' @return a data frame
#' @details // TO DO. --> explain columns 'signif' and 'type' in the returned data frame.
#' @examples
#' \dontrun{
#' // TO DO
#' }
categorize_diff_genes <-  function(dsr_tibble,
                                   log2fc_cutoff =  1,
                                   pval_cutoff = 0.05,
                                   padj_cutoff = 0.01,
                                   add_column_regul = TRUE,
                                   regul_based_upon = 1
                                   ){

  # convert DESeqResult object to tibble

  if(any(is(dsr_tibble) %in% "DESeqResults")) {
    dsr_tibble %<>% dsr_to_tibble()
  }

  # add column 'type' having one of the four values

  # 1. NS -> Not significant
  # 2. p-value -> only pvalue significant but gene is not diff expressed
  # 3. log2FC -> only diff expressed but not pvalue significant
  # 4. p-value&log2FC -> diff expressed and pvalue signficant

  pval_cutoff <- pval_cutoff
  log2fc_cutoff <- log2fc_cutoff

  dsr_tibble %<>% dplyr::mutate(signif = dplyr::case_when(dplyr::between(log2FoldChange , -log2fc_cutoff, log2fc_cutoff ) &  pvalue >= pval_cutoff ~ "NS",
                                                        dplyr::between(log2FoldChange , -log2fc_cutoff, log2fc_cutoff  ) &  pvalue <= pval_cutoff ~ "p-value",
                                                        (log2FoldChange <= -log2fc_cutoff | log2FoldChange >= log2fc_cutoff) & pvalue >= pval_cutoff ~ "log2FC",
                                                        (log2FoldChange <= -log2fc_cutoff | log2FoldChange >= log2fc_cutoff) & pvalue <= log2fc_cutoff ~ "p-value&log2FC",
                                                        TRUE ~ "other"))

  # add column 'regul' having one of the three values

  if(add_column_regul){

    # categorization into 'Up', 'Down' and 'Other' is done based on value of user selected value of regul_based_upon

    # if regul_based_upon == 1, filtering is done based on pvalue cutoff
    #   Up -> log2fc >= log2fc_cutoff & pvalue <= pvalue_cutoff
    #   Down -> log2fc  <= (-1) * log2fc_cutoff & pvalue <= pvalue_cutoff
    #   Other -> remaining genes

    if(regul_based_upon == 1){
      dsr_tibble %<>%
        dplyr::mutate(regul = dplyr::case_when(log2FoldChange >= log2fc_cutoff & pvalue <= pval_cutoff ~ "Up",
                                               log2FoldChange <= -1 * log2fc_cutoff & pvalue <= pval_cutoff ~ "Down",
                                               TRUE ~ "other"))

    }
    # if regul_based_upon == 2

    #   Up -> log2fc >= log2fc_cutoff & padj <= padj_cutoff
    #   Down -> log2fc  <= (-1) * log2fc_cutoff & padj <= padj_cutoff
    #   Other -> remaining genes

    else if(regul_based_upon == 2){
      dsr_tibble %<>%
        dplyr::mutate(regul = dplyr::case_when(log2FoldChange >= log2fc_cutoff & padj <= padj_cutoff ~ "Up",
                                               log2FoldChange <= -1 * log2fc_cutoff & padj <= padj_cutoff ~ "Down",
                                               TRUE ~ "other"))

    }

    # if regul_based_upon == 3

    #   Up -> log2fc >= log2fc_cutoff & pvalue <= pvalue_cutoff & padj <= padj_cutoff
    #   Down -> log2fc  <= (-1) * log2fc_cutoff pvalue <= pvalue_cutoff & padj <= padj_cutoff
    #   Other -> remaining genes

    else if(regul_based_upon == 3){
      dsr_tibble %<>%
        dplyr::mutate(regul = dplyr::case_when(log2FoldChange >= log2fc_cutoff &
                                                 pvalue <= pval_cutoff &
                                                 padj <= padj_cutoff   ~ "Up",
                                               log2FoldChange <= -1 * log2fc_cutoff &
                                                 pvalue <= pval_cutoff &
                                                 padj <= padj_cutoff ~ "Down",
                                               TRUE ~ "other"))

    } else {
      stop ("value of 'regul_based_upon' must be one of the 1,2, or 3.")
    }
  }

}



#' Prepare a fold change matrix
#' @description This function returns a dataframe having first column gene names and subsequent columns are
#' fold change values for the comparisons passed through `sample_comparisons`.
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparisons a character vector denoting sample comparisons for which fold change values to be derived.
#' @param genes a character vector denoting gene names for which fold change values to be derived.
#'
#' @return a tbl.
#' @export
#'
#' @examples
#' \dontrun{
#' // TO DO
#' }
#'
get_fold_change_matrix <- function(x , sample_comparisons , genes){

  # validate x

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))

  # validate sample_comparisons

  stopifnot("sample_comparisons must be a character vector." = is.character(sample_comparisons))

  # validate genes

  stopifnot("genes must be a character vector." = is.character(genes))

  # check if sample_comparisons present in x

  if(!any (x$comp %in% sample_comparisons)) {
    stop("None of values from 'sample_comparisons' present in the x.")
  }

  # get gene_id column name

  ## NOTE: There should be a way to get all user input parameters from the 'parcutils' object.

  gene_id_column <- x$dsr_tibble_deg[[1]] %>% colnames() %>% .[1]

  # subset data
  res <- x %>%
    dplyr::filter(comp %in% sample_comparisons ) %>%
    dplyr::pull(dsr_tibble_deg) %>%
    purrr::map(~ ..1 %>% dplyr::filter(!!rlang::sym(gene_id_column) %in% genes) %>% dplyr::select(1,3)) %>%
    dplyr::bind_rows(.id  = "comparisons") %>%
    tidyr::pivot_wider(id_cols = !!rlang::sym(gene_id_column) , names_from = comparisons, values_from = log2FoldChange)


  if(nrow(res) == 0){
    warning("None of the queried 'genes' are found in the 'x'.")
  }

  return(res)
}



#' Prepare a matrix of normalised gene expression values
#' @description This function returns a dataframe having first column gene names and subsequent columns are
#' normalised gene expression values for the comparisons passed through sample_comparisons.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting samples for which normalised gene expression values to be derived.
#' @param genes a character vector denoting gene names for which normalised gene expression values to be derived.
#' @param summarise_replicates logical, default FALSE, indicating whether gene expression values summarised by mean or median between replicates.
#' @param summarise_method a character string either "mean" or "median" by which normalised gene expression values will be summarised between replicates.
#'
#' @return a tbl.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' // TO DO
#' }
get_normalised_expression_matrix <- function(x , samples, genes, summarise_replicates = FALSE, summarise_method = "median" ){

  # validate x

  validata_parcutils_obj(x)

  # validate genes

  stopifnot("genes must be a character vector." = is.character(genes))

  # validate summarise_replicates

  stopifnot("summarise_replicates must be a logical." = is.logical(summarise_replicates))

  # validate samples

  stopifnot("samples must be a character vector." = is.character(samples))

  # validate summarise_method

  match.arg(summarise_method , choices = c("mean" ,"median"))

  gene_id_column <- x$norm_counts[[1]][[1]] %>% colnames() %>% .[1]


  ## get all samples gene expression matrix.
  all_expr_mats <-  get_all_named_expression_matrix(x)

  # merge df and make them long format
  expr_mat_long <- all_expr_mats[samples] %>%
    purrr::map_df(~ ..1 %>% tidyr::pivot_longer(cols = -1 , values_to = "vals" , names_to = "replicate") , .id = "sample") %>%
    dplyr::select(-1,dplyr::everything(),1) %>%
    dplyr::distinct()

  column_gene_id <- expr_mat_long %>% colnames()%>%.[1]

  # summarise replicates

  if(summarise_replicates){
    expr_mat_long <- expr_mat_long %>%
      dplyr::group_by(sample , !!rlang::sym(column_gene_id))  %>%
      dplyr::summarise_at("vals", summarise_method) %>%
      dplyr::ungroup()

    #  make it wide
    expr_mat_wide <- expr_mat_long %>%
      tidyr::pivot_wider(id_cols = !!rlang::sym(column_gene_id), values_from = "vals", names_from = "sample")

    # keep original order of columns
    expr_mat_wide <- expr_mat_wide %>%
      dplyr::select(!!rlang::sym(column_gene_id),dplyr::all_of(samples))

  } else{
    #  make it wide
    expr_mat_wide <- expr_mat_long %>%
      tidyr::pivot_wider(id_cols = !!rlang::sym(column_gene_id), values_from = "vals", names_from = "replicate" )

  }

  # filter by user supplied genes

  expr_mat_wide <- filter_df_by_genes(df = expr_mat_wide, genes = genes)

  return(expr_mat_wide)


}



#' Get genes from based on their differential regulation (up, down, both, other and all)
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparison a character string denoting a sample comparison for which genes to be obtained.
#' @param regulation a character string, default \code{both}. Values can be one of the \code{up}, \code{down}, \code{both}, \code{other}, \code{all}.
#'  + up : returns all up regulated genes.
#'  + down : returns all down regulated genes.
#'  + both : returns all up and down regulated genes.
#'  + other : returns genes other than up and down regulated genes.
#'  + all : returns all genes.
#' @return a named vector.
#' @export
#'
#' @examples
#' \dontrun{
#' // TO DO
#' }
#'
get_genes_by_regulation <-  function(x, sample_comparison , regulation = "both"  ) {

  # validate x.
  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))

  # validate sample comparison.
  stopifnot("sample_comparison must be a character string" = is.character(sample_comparison) & length(sample_comparison) ==1)

  # validate regulation.
  match.arg(regulation , choices = c("up","down" , "both", "other" ,"all"), several.ok = F)

  # split by column regul.
  # NOTE: do not group by column name instead use index. The reason is because if the column name of column 'regul' change in future it will break this code.

  genes_by_comp <- x$dsr_tibble_deg[[sample_comparison]] %>%
    dplyr::select(1, dplyr::last_col())

  genes_by_comp <- genes_by_comp %>%
    dplyr::group_by_at(2)

  grp_keys <- genes_by_comp %>%
    dplyr::group_keys() %>%
    pull(1)

  genes_by_comp <- genes_by_comp %>%
    dplyr::group_split()
  names(genes_by_comp) <- grp_keys

  # Convert names and regulation to lower case before comparison
  names(genes_by_comp)  <- tolower(names(genes_by_comp))
  regulation = tolower(regulation)

  if(regulation == "both"){
    rslt <- genes_by_comp[c("up","down")]  %>% as.list() %>% dplyr::bind_rows()
  } else if(regulation == "all"){
    rslt <- genes_by_comp  %>% as.list() %>% dplyr::bind_rows()
  } else{
    rslt <- genes_by_comp[[regulation]]
  }

  rslt <- rslt[[1]] %>%
    purrr::set_names(rslt[[2]] %>% tolower())

  return(rslt)

}







#' Generate upset plots for DEG between comparisons.
#' @description  Given a set of DEG comparisons, the functions returns [UpSetR::upset()] plots for up and down genes for any 2 comparisons.
#' For each upset plot generated function also returns interaction in form of tibble.
#' @param x x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparisons a character vector denoting  sample comparisons between upset plot to be generated.
#' @param color_up a character string denoting a valid color code for bars in upset plot for up regulated genes.
#' @param color_down a character string denoting a valid color code for bars in upset plot for down regulated genes.
#'
#' @return an object of named list where each element is a list of two - 1)  an upset plots  [UpSetR::upset()] and their intersections in form of tibble.
#' @export
#' @importFrom  UpSetR upset fromList
#' @importFrom  purrr map set_names cross map_chr
#' @importFrom  glue glue
#' @importFrom  magrittr %>% %<>%
#' @examples
#' \dontrun{
#'
#' // TO DO.
#' }
#'
plot_deg_upsets <- function(x, sample_comparisons, color_up = "#b30000", color_down = "#006d2c"){

  # x <- dd
  # sample_comparisons = c("shRNA1_VS_shCTRL", "INFy_shRNA1_VS_INFy_shCTRL" , "INFy_shCTRL_VS_shCTRL")

  # validate x

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))


  # validate sample_comparisons

  stopifnot("sample_comparisons must be a character vector of length 2." = is.character(sample_comparisons) & length(sample_comparisons) > 1)

  # all sample comparisons must present in x

  if(!all(sample_comparisons %in% x$comp) ){
    not_present <- sample_comparisons[!(sample_comparisons %in% x$comp)]
    stop(glue::glue("Value {not_present} is not present in the x. Check x$comp to access all choices."))
  }

  # validate colors

  stopifnot("color_up and color_down must be a character string denoting a valid color code." = is.character(color_up) & is.character(color_up) & length(color_up)==1, length(color_down)==1)

  # generate all combinations of 2 from sample_comparisons.
  sample_comparisons_all_comb <- combn(sample_comparisons, 2)  %>%
    as.data.frame() %>%
    purrr::map( ~..1)


  ## set names
  sample_comparisons_all_comb <- purrr::set_names(sample_comparisons_all_comb , ~purrr::map_chr(sample_comparisons_all_comb , ~paste0(..1, collapse = "_AND_")))


  all_upsets <- purrr::map(sample_comparisons_all_comb, ~ piarwise_upset(x = x,
                                                                         sample_comparison = ..1,
                                                                         color_up = color_up ,
                                                                         color_down = color_down) )

  return(all_upsets)
}




##
#' Generate upset plots for DEG between a comparison.
#' @description This function called from [parcutils::plot_deg_upsets()].
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparison a character vector of length 2 denoting  sample comparisons between upset plot to be generated.
#' @param color_up a character string denoting a valid color code for bars in upset plot for up regulated genes.
#' @param color_down a character string denoting a valid color code for bars in upset plot for down regulated genes.
#' @param ... other arguments to be passed to [UpSetR::upset()].
#'
#' @return an object of named list where each element is a list of two - 1)  an upset plots  [UpSetR::upset()] and their intersections in form of tibble.
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#'
#' // TO DO.
#' }
#'
piarwise_upset <- function(x, sample_comparison , color_up = "#b30000", color_down = "#006d2c",... ){

  # validate sample_comparison

  stopifnot("sample_comparison must be a character vector of length 2." = is.character(sample_comparison) & length(sample_comparison) == 2)

  yy <- purrr::map(sample_comparison %>%
                     purrr::set_names(sample_comparison), ~parcutils::get_genes_by_regulation(x = x, sample_comparison = ..1,regulation = "both") ) %>%
    purrr::map(~split(..1, names(..1)))

  input_list_for_upset <- purrr::map(names(yy), ~ purrr::set_names(yy[[..1]], glue::glue("{..1}_{yy[[..1]] %>% names %>% tolower}"))) %>% purrr::flatten()


  # order list by up and down genes

  lst_order <- list(sample_comparison , c("up" ,"down")) %>%
    purrr::cross() %>%
    purrr::map_chr( paste0 , collapse ="_")

  input_list_for_upset <- input_list_for_upset[lst_order]

  # up set colors
  upset_colors <- c("up" = color_up, #MetBrewer::met.brewer(name = "Austria")[1],
                    "down" = color_down) #MetBrewer::met.brewer(name = "Austria")[3])



  # plot upset
  us_plot <- UpSetR::upset(UpSetR::fromList(input_list_for_upset) ,
                           order.by =  c("degree"),
                           #group.by = "sets",
                           nsets = length(input_list_for_upset),
                           nintersects = NA,
                           sets = names(input_list_for_upset) %>% rev(),
                           keep.order=T,
                           text.scale=1.5,
                           mb.ratio = c(0.6, 0.4),
                           sets.x.label = "No. of genes",
                           sets.bar.color = rep(upset_colors %>% rev(),each = 2),...)

  # get upset data.
  upset_intersects  <- parcutils::get_upset_intersects(input_list_for_upset, upset_plot = us_plot)

  ret <- list(upset_plot = us_plot,upset_intersects = upset_intersects)


  return(ret)


}




#' Generate a heatmap of normalised gene expression values, z-score or log2 fold-change values.
#' @description Heatmap is a common tool to show gene expression pattern across samples in RNA-seq experiment.
#' While doing data exploration, it is a common practice to generate several heatmaps to identify interesting
#' patterns of gene expression across sample. Most heatmap generating tools require data in a tabular form.
#' However, to prepare such a tabular form requires significant data wrangling such as sub-setting
#' genes (rows) and samples (columns) specifically when RNA-seq studies studies consisting of several samples and
#' replicates.
#'
#' This function cut downs several steps of data wrangling to create a heatmap of gene expression, z-score or
#' log2 fold-change . Given an object of 'parcutils', names of samples or sample comprisons, genes to show in heatmap
#' and several other arguments it creates a heatmap.The output of the function is an output of the function [ComplexHeatmap::Heatmap()]
#' which then can be used for the other functions of the [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/) package.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting sample names to show in the heatmap.
#' @param genes a character vector denoting genes to show in the heatmap.
#' @param repair_genes  logical, default FALSE, indicating whether to repair gene names or not. See details.
#' @param convert_log2 logical, default FALSE, indicating whether to convert gene expression values in log2 or not.
#' @param color_default logical, default TRUE, indicating whether to use default heatmap colors.
#' @param col an output of [circlize::colorRamp2()], default NULL.
#' @param convert_zscore logical, default TRUE, indicating whether to convert gene expression values in to the z-score or not. see details.
#' @param summarise_replicates logical, default TRUE, indicating whether to summarise values for each gene across replicates.
#' @param summarise_method a character string, either mean or median.
#' @param show_row_names logical, default FALSE, indicating whether to show row names in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param cluster_rows logical, default TRUE, indicating whether to cluster rows in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param show_row_dend logical, default TRUE, indicating whether to show dendrogram in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param row_names_font_size a numeric value, default 10, indicating size of row names in the heatmap.
#' @param show_column_names logical, default TRUE, indicating whether to show column names in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param cluster_columns logical, default TRUE, indicating whether to cluster columns in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param show_heatmap_legend logical, default TRUE, indicating whether to show heatmap legend or not.
#' @param ... Other parameters to be passes to [ComplexHeatmap::Heatmap()].
#'
#' @return an output of the function [ComplexHeatmap::Heatmap()].
#' @export
#' @details
#' The function `get_gene_expression_heatmap` is to create a heatmap either for normalised gene expression or z-score values
#' while the function `get_fold_change_heatmap` is to create a heatmap for log2 fold-changes values.
#'
#' + \code{repair_genes} :  Internally gene names are stored as a "gene_id:gene_symbol" format. For example, "ENSG00000187634:SAMD11".
#' When \code{repair_genes} is set to \code{TRUE} the string corresponding to gene_id followed by ":" will be removed. This is useful when gene names
#' to be revealed in the heatmap.
#'
#' + \code{convert_zscore} :  When set to \code{TRUE} values for each gene is converted into z-score. z-score is calculated by baser r function [base::scale()]
#' with all default parameters.
#'
#' @examples
#' \dontrun{
#' // TO DO.
#' }
#'
get_gene_expression_heatmap <- function(x,
                                        samples,
                                        genes,
                                        repair_genes = FALSE,
                                        convert_log2 = FALSE,
                                        color_default = TRUE,
                                        col = NULL,
                                        convert_zscore = TRUE,
                                        summarise_replicates = TRUE,
                                        summarise_method = "median",

                                        show_row_names = FALSE,
                                        cluster_rows = TRUE,
                                        show_row_dend = TRUE,
                                        row_names_font_size = 10,

                                        show_column_names = TRUE,
                                        cluster_columns = TRUE,

                                        show_heatmap_legend = TRUE,

                                        ...){


  ## argument row_km from ComplexHeatmap::Heatmap is not supported.
  arg_dots <- list(...)
  if("row_km" %in% names(arg_dots)){
    stop("Due to inconsistant output of ComplexHeatmap::row_order() between arguments `row_split` and `row_km`, `row_km` is currently not supported. Use `row_split` instead.")
  }


  ## validate x

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))

  # validate samples

  stopifnot("samples must be a character vector." = is.character(samples))

  # validate genes

  stopifnot("genes must be a character vector." = is.character(genes))

  # validate repair_genes

  stopifnot("repair_genes must be a logical." = is.logical(repair_genes))


  # validate color default

  stopifnot("color_default must be a logical." = is.logical(color_default))

  if(!color_default){
    if(is.null(col)){
      stop("col must be the output of circlize::colorRamp2() when  color_default is FALSE.")
    }
  }

  # validate summarise_replicates

  stopifnot("summarise_replicates must be a logical." = is.logical(summarise_replicates))

  # validate summarise_method

  match.arg(summarise_method , choices = c("mean" ,"median"))

  # validate convert_zscore

  stopifnot("convert_zscore must be a logical." = is.logical(convert_zscore))

  ## prepare matrix for heatmap
  all_expr_mat <-  get_all_named_expression_matrix(x)

  # merge df and make them long format
  expr_mat_long <- all_expr_mat[samples] %>%
    purrr::map_df(~ ..1 %>% tidyr::pivot_longer(cols = -1 , values_to = "vals" , names_to = "replicate") , .id = "sample") %>%
    dplyr::select(-1,dplyr::everything(),1) %>%
    dplyr::distinct()

  column_gene_id <- expr_mat_long %>% colnames()%>%.[1]

  # summarise replicates

  if(summarise_replicates){
    expr_mat_long <- expr_mat_long %>%
      dplyr::group_by(sample , !!rlang::sym(column_gene_id))  %>%
      dplyr::summarise_at("vals", summarise_method) %>%
      dplyr::ungroup()

    #  make it wide
    expr_mat_wide <- expr_mat_long %>%
      tidyr::pivot_wider(id_cols = !!rlang::sym(column_gene_id), values_from = "vals", names_from = "sample")

    # keep original order of columns
    expr_mat_wide <- expr_mat_wide %>%
      dplyr::select(!!rlang::sym(column_gene_id),dplyr::all_of(samples))

  } else{
    #  make it wide
    expr_mat_wide <- expr_mat_long %>%
      tidyr::pivot_wider(id_cols = !!rlang::sym(column_gene_id), values_from = "vals", names_from = "replicate" )

  }

  # convert log
  if(convert_log2){
    expr_mat_wide <- expr_mat_wide %>% TidyWrappers::tbl_convert_log2(frac = 0.1)
  }

  # convert z score

  if(convert_zscore){
    expr_mat_wide <- expr_mat_wide %>% TidyWrappers::tbl_convert_row_zscore()
  }

  # filter by user supplied genes

  expr_mat_wide <- filter_df_by_genes(df = expr_mat_wide, genes = genes)

  # remove row if all are NA.
  value_na <- expr_mat_wide  %>%  TidyWrappers::tbl_keep_rows_NA_any() %>%
    dplyr::pull(!!rlang::sym(column_gene_id))
  if(length(value_na) > 0){
    value_na <- paste0(value_na, collapse = ",")
    warning(glue::glue("Genes having value NA - {value_na} are removed from heatmap.",))
    expr_mat_wide <- expr_mat_wide %>% TidyWrappers::tbl_remove_rows_NA_any()
  }

  #fix colors
  if(color_default){
    col = fix_hm_colors(expr_mat_wide)
  } else {
    col = col
  }

  # repair genes

  if(repair_genes){
    expr_mat_wide <- expr_mat_wide %>% dplyr::mutate(!!column_gene_id := !!rlang::sym(column_gene_id) %>% stringr::str_replace(".*:" , ""))
  }

  #  convert into matrix
  expr_mat_wide  %<>%
    as.data.frame() %>% tibble::column_to_rownames(column_gene_id) %>% as.matrix()

  # generate heatmap

  hm <- ComplexHeatmap::Heatmap(expr_mat_wide,col = col,
                                show_row_names = show_row_names,
                                cluster_rows = cluster_rows,
                                show_row_dend = show_row_dend,
                                row_names_gp = grid::gpar(fontsize = row_names_font_size),
                                show_column_names = show_column_names,
                                cluster_columns = cluster_columns,
                                show_heatmap_legend = show_heatmap_legend, ...)

  return(hm)

}


#' @rdname get_gene_expression_heatmap
#'
#' @param sample_comparisons a character vector denoting sample comparisons for which heatmap of log2 fold-change to be plotted.
#' @export
get_fold_change_heatmap <-  function(x,
                                     sample_comparisons,
                                     genes,
                                     repair_genes = FALSE,
                                     color_default = TRUE,
                                     col = NULL,

                                     show_row_names = FALSE,
                                     cluster_rows = TRUE,
                                     show_row_dend = TRUE,
                                     row_names_font_size = 10,

                                     show_column_names = TRUE,
                                     cluster_columns = TRUE,

                                     show_heatmap_legend = TRUE, ...){

  arg_dots <- list(...)
  if("row_km" %in% names(arg_dots)){
    stop("Due to inconsistant output of ComplexHeatmap::row_order() between arguments `row_split` and `row_km`, `row_km` is currently not supported. Use `row_split` instead.")
  }


  ## validate x

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))

  # validate samples

  stopifnot("sample_comparisons must be a character vector." = is.character(sample_comparisons))

  # validate genes

  stopifnot("genes must be a character vector." = is.character(genes))

  # validate repair_genes

  stopifnot("repair_genes must be a logical." = is.logical(repair_genes))


  # validate color default

  stopifnot("color_default must be a logical." = is.logical(color_default))

  if(!color_default){
    if(is.null(col)){
      stop("col must be the output of circlize::colorRamp2() when  color_default is FALSE.")
    }
  }


  fc_df <-  parcutils::get_fold_change_matrix(x = x,
                                              sample_comparisons = sample_comparisons,
                                              genes = genes)

  gene_id_column <- fc_df %>% colnames() %>% .[1]

  # check if any of the given gene not found

  genes_not_found <- genes[! genes %in% fc_df[[1]] ]

  if(length(genes_not_found) > 0 ){
    cli::cli_alert_warning(paste0("Below genes are not found in x.",  paste0(genes_not_found , collapse = ",")))
  }

  # keep genes and samples in the same order

  fc_df <- fc_df %>%
    dplyr::select(!!rlang::sym(gene_id_column),dplyr::all_of(sample_comparisons)) %>%
    dplyr::slice(match(genes, fc_df[[1]]))


  # fix colors
  if(color_default){
    col = fix_hm_colors(fc_df)
  } else {
    col = col
  }

  # repair genes

  if(repair_genes){
    fc_df <- fc_df %>%
      dplyr::mutate(!!gene_id_column := !!rlang::sym(gene_id_column) %>%
                      stringr::str_replace(".*:" , ""))
  }

  # convert df to matrix
  fc_mat <- fc_df %>% as.data.frame() %>% tibble::column_to_rownames(gene_id_column) %>% as.matrix()

  # generate heatmap

  hm <- ComplexHeatmap::Heatmap(fc_mat,
                                col = col,
                                show_row_names = show_row_names,
                                cluster_rows = cluster_rows,
                                show_row_dend = show_row_dend,
                                row_names_gp = grid::gpar(fontsize = row_names_font_size),
                                show_column_names = show_column_names,
                                cluster_columns = cluster_columns,
                                show_heatmap_legend = show_heatmap_legend, ...)

  return(hm)
}



#' Filter expression matrix (dataframe) by genes(first column).
#'
#' @param df a gene expression data frame
#' @param genes acharacter vector of genes to be filtred from df.
#'
#' @return a filtered dataframe.
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#'
#' # TO DO
#' }
#'
filter_df_by_genes <- function(df, genes){

  column_gene_id <- df %>% colnames() %>%.[1]

  ## check which genes are not present in the df
  df_genes <- df %>%
    dplyr::pull(!!rlang::sym(column_gene_id))

  genes_not_present <- dplyr::setdiff(genes, df_genes)
  genes_present <- genes[genes %in% df_genes]

  if(length(genes_not_present) == length(genes) ){

    stop("No genes found.")

  } else if(length(genes_present) < length(genes)){

    cli::cli_alert_warning(text = paste0("Genes below are not present into x\n" ,
                                         paste0(genes_not_present , collapse = ",")))

  }

  # filter available genes
  df <- df %>%
    dplyr::slice(match(genes_present, !!rlang::sym(column_gene_id)))

}


#' Get named list of gene expression matrix for all samples in x.
#'
#' @param x an object of class parcutils.
#'
#' @return named list of gene expression matrix.
#' @export
#' @keywords internal
get_all_named_expression_matrix <- function(x){

  validata_parcutils_obj(x)

  all_sample_names <- x$norm_counts %>%
    purrr::flatten() %>% names()

  unique_index <- purrr::map_dbl(all_sample_names %>% unique() , ~ which(..1 == all_sample_names)[1])

  ret <- x$norm_counts %>% purrr::flatten() %>% magrittr::extract(unique_index)
  return(ret)

}


#' Validate if an object is of class 'parcutils'.
#'
#' @param x
#'
#' @return
#' @export
#' @keywords internal
validata_parcutils_obj <- function(x){

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))
}



#' Given a matrix generate color scale for heatmap.
#'
#' @param hm_matrix
#'
#' @return color for [ComplexHeatmap::Heatmap()].
#' @export
#' @keywords internal
fix_hm_colors <- function(hm_matrix){

  stopifnot(is.data.frame(hm_matrix))

  col_choices <- c("min" = MetBrewer::met.brewer(name = "Austria")[2],
                   "max" = MetBrewer::met.brewer(name = "Austria")[1])

  max_val <- hm_matrix %>% dplyr::select_if(is.numeric) %>% max() %>% ceiling()
  min_val <- hm_matrix %>% dplyr::select_if(is.numeric) %>% min() %>% floor()
  mid_val <- (max_val + min_val) / 2

  col = circlize::colorRamp2(c(min_val , mid_val, max_val), c(col_choices["min"], "white", col_choices["max"]))


}



## given a heatmap,  x and destination file  save expression values and / or zscore values.

#' Get data from a heatmap in the same order.
#'
#' @param h an object of class [ComplexHeatmap::Heatmap()]
#'
#' @return a tbl.
#' @export
#'
#' @examples
#' \dontrun{
#' # // TO DO
#' }
#'
get_heatmap_data <- function(h){


  # validate h
  if(!is(h , "Heatmap")){
    stop("h must be the output of ComplexHeatmap::Heatmap().")
  }

  h <- ComplexHeatmap::draw(h)

  # get row order.
  row_ord <- ComplexHeatmap::row_order(h)

  # if not list, make it.

  if(!is(row_ord, "list")){
    row_ord <- list(row_ord)
  }

  # set names
  row_ord <- row_ord %>% purrr::set_names(nm = glue::glue("cluster_{1:length(row_ord)}"))

  ## make tibble
  row_ord <- tibble::tibble(clust = names(row_ord) , ord = row_ord) %>% tidyr::unnest(cols = (ord))

  mat <- h@ht_list[[1]]@matrix %>% as.data.frame() %>% tibble::rownames_to_column()  %>%
    tibble::as_tibble()

  ## order matrix by hm
  mat_ordered <- mat %>%
    dplyr::slice(row_ord$ord)

  ## add clusters
  mat_ordered <- mat_ordered %>%
    dplyr::mutate(clusters = row_ord$clust)

  return(mat_ordered)
}











