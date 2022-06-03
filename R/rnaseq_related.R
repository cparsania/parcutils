#' Run differential expression (DE) analysis for several comparisosns at a time. Internally it uses
#' [DESeq2::DESeq()] to perform DE analysis.
#'
#' @description
#'
#' DESeq2 is a popular method to perform DE analysis for RNA-seq data. Pre and post DESeq run, however,
#' involves several data wrangling steps. For example, prior to run DESeq genes may need to be filtered out based on
#' number of reads mapped to genes across replicates. Similarly,  post DESeq run user needs to set cutoffs (log2fc, pvalue and padj)
#' to define `up` and `down` regulated genes. For large RNA-seq experiments involving several DE comparison such as sub-setting
#' and different cutoffs creates messy and less readable code. This function helps to make things bit tidy and increase
#' the code readability. Besides that,  output of this function can subsequently used for several other functions of this package.
#'
#' @param counts a character string of the path to a count file or an object of dataframe having raw counts for each gene.
#' See details below to know more about format of the count file and count dataframe.
#' @param column_geneid a character string denoting a column of geneid in \code{counts}
#' @param sample_info a character string denoting a name of sample information file or a dataframe.
#' A file or a dataframe both must have at least two columns *WITHOUT* column names. First column denotes to samples names
#' and second column denotes group name for each sample in first column. For e.g.
#'
#'  |       |  |
#'  | ----------- | ----------- |
#'  | control_rep_1      | Control       |
#'  | control_rep_2      | Control       |
#'  | control_rep_3      | Control       |
#'  | treat_1_rep_1     | Treatment_1       |
#'  | treat_1_rep_2      | Treatment_1       |
#'  | treat_1_rep_3      | Treatment_1       |
#'  | treat_2_rep_1      | Treatment_2      |
#'  | treat_2_rep_2      | Treatment_2       |
#'  | treat_2_rep_3      | Treatment_2       |
#'
#' @param group_numerator a character vector denoting sample groups to use in numerator to calculate fold change.
#' @param group_denominator a character vector denoting sample groups to use in denominator to calculate fold change.
#' @param delim a character denoting deliminator for `count` file. Only valid if `count` is a file path.
#' @param comment_char a character denoting comments line in count file. Only valid if `count` is a file path.
#' @param min_counts a numeric value, default 10,  denoting minimum counts for a gene to be used to consider a
#' gene for differential expression analysis.
#' @param min_replicates a numeric value, default 1, denoting minimum samples within a sample group must have `minimum_counts`.
#' Value provided must not be higher than number of samples in each group.
#' For example for given values `min_replicates = 2` and  `minimum_counts = 10`
#' the genes which have minimum counts 10 in at least 2 samples within a groups will be accounted for DEG. Rest will be filtered out prior to run DESeq.
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
#' @param cutoff_lfc minimal threshold for log2fold change, default 1 (2 fold).
#' @param cutoff_pval minimal threshold for pvalue, default 0.05. P-value threshold will be applied only when
#'  `regul_based_upon` is either 1 or 3.
#' @param cutoff_padj minimal threshold for Padj, default 0.01. Padj threshold will be applied only when
#'  `regul_based_upon` is either 2 or 3.
#' @param print_rows_all_zero logical, default FALSE, denoting whether to print genes with value 0 in all columns.
#'
#' @details : For the argument `count` user can either provide a character string denoting a file or an object of dataframe.
#'  In each case required format is explained below.
#' + `Count file`: Count file is a table of row read counts usually derived from .bam file for different genomic features
#'  (e.g. genes, transcripts etc.). Data in the count file must be in a tabular format with a valid column deliminator
#'  (e.g., tab, comma etc.). First row and first column will be considered as column names and row names respectively.
#'  Values for column names and row names are usually character string or combination of character, numbers and
#'  special characters such as `_`, or  `.`. Both row names and column names must have unique values.
#' + `Count dataframe`: Count data in a dataframe format having same requirement of row names and column names explained for `count file`.
#' @return Return object is a dataframe (tibble) having each row denoting a unique differential comparison. There are total 8 columns as explained below.
#' + `de_comparisons` : It stores the name of differential comparison for each row.
#' + `numerator` : It stores name of samples which were used as numerator for the differential comparison in each row.
#' + `denominator` : It stores name of samples which were used as denominator for the differential comparison in each row.
#' + `norm_counts` : This is a `named-list` column. Each row in this column is a list of two containing normalised
#' genes expression values in a dataframe for the samples - numerator and denominator. The first column of the dataframe is `gene_id`
#' and subsequent columns are gene expression values in replicates of corresponding samples. This normalised gene expression values are
#' obtained using  `counts` slot of a DESeqDataSet object. e.g.: `counts(dds,normalized=TRUE)`
#' + `dsr` : This is a `named-list` column stores an object of class [DESeq2::DESeqResults()] for the
#' differential comparison in each row.
#' + `dsr_tibble`:  This is a `named-list`column stores and an output of [DESeq2::DESeqResults()] in the dataframe format for the
#' differential comparison in each row.
#' + `dsr_tibble_deg`: The data in this column is same as in the column `dsr_tibble` except it contains two extra columns `signif` and `regul`.
#' Values in the `signif` specifies statistical and fold change significance of the gene while values in the `regul` denotes whether the gene is `up` or `down`
#' regulated.
#' + `deg_summmary` : This is a `named-list` column. Each element of the list is a dataframe summarizing number of differential expressed gene for the differential
#' comparison for each row.
#' @export
#' @importFrom tidyselect everything
#' @examples
#'
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t",show_col_types = FALSE)
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          cutoff_lfc = 1 ,
#'                          cutoff_pval = 0.05,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#'
#' res
#'
#' ## all comparisons
#'
#' print(res$de_comparisons)
#'
#' ## DESEq result object(s)
#'
#' print(res$dsr)
#'
#' ## DESEq result data frame
#'
#' print(res$dsr_tibble)
#'
#' ## DESEq result data frame  DEG assigned, look at the columns 'signif' and 'regul'
#'
#' print(res$dsr_tibble_deg)
#'
#' ## DEG summary
#'
#' print(res$deg_summmary)
run_deseq_analysis <- function(
  counts,
  column_geneid,
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
  print_rows_all_zero = FALSE
){

  # defaults
  # counts <- matrix(rnbinom(n=1000, mu=100, size=1/0.5), ncol=10) %>% as.data.frame() %>% tibble::as_tibble()
  #
  # colnames(counts) <- c(paste("c" , c(1:5), sep = ""),c(paste ("d" , 1:5, sep = "")))
  #
  # counts %<>%  dplyr::mutate("Geneid" = stringi::stri_rand_strings(n = 100, length = 5))  %>% dplyr::relocate("Geneid")
  # sample_info <- tibble::tibble(samples = colnames(counts)[-1] , sample_groups = factor(rep(c("c","d"), each=5)))
  # column_geneid = "Geneid"
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
      cli::cli_abort("{names(is_numeric_and_scaler[i])} must be of length one and of class {.cls double}.")
  }

  if(!rlang::is_scalar_double(min_counts)){
    cli::cli_abort(cli::cli_text("{.arg min_counts} must be of length one of class {.cls double}."))
  }

  # value for an argument column_geneid must be a character string

  if(!is.character(column_geneid) && length(column_geneid) != 1){
    cli::cli_abort("Value for an argument {.arg geneid_coulmn} must be a {.cls character} string")
  }

  # quote column names
  column_geneid_quo <- rlang::enquo(column_geneid)

  ## prepare count matrix
  count_data <- .get_count_data(counts)

  # prepare sample_info
  sample_info_data <- .get_sample_information(sample_info)

  # assign column names to sample_info_data
  colnames(sample_info_data) <- sample_info_colnames

  # check if all the samples provided in sample_info are present in count_data file
  count_data_col_names <-  colnames(count_data)
  sample_info_samples <- sample_info_data %>% dplyr::pull(1)
  sample_info_samples_quo <- rlang::enquos(sample_info_samples)

  if(!all(sample_info_samples %in% count_data_col_names)){
    index_not_present <- which(!(sample_info_samples %in% count_data_col_names))
    value_not_present <- sample_info_samples[index_not_present] %>%
      paste0(collapse = "','")
    cli::cli_abort("Values - {.emph {cli::col_red({value_not_present})}} -
                   provided in {.arg sample_info} (first column) must
                   present as column{?s} in {.arg counts}.")
  }

  # validate 'group_numerator' and 'group_denominator'

  if(!inherits(group_numerator , what = "character")) {
    cli::cli_abort(" {.arg group_numerator} must be a {.cls character} vector.")
  }

  if(!inherits(group_denominator , what = "character")) {
    cli::cli_abort("{.arg group_denominator} must be a {.cls character} vector.")
  }

  # values mentioned for 'group_numerator' and 'group_denominator' must present in sample_info
  sample_info_groups <- sample_info_data %>%
    dplyr::pull(sample_info_colnames[2]) %>%
    unique()

  if(!all(group_numerator %in% sample_info_groups)){
    cli::cli_abort("all values from {.arg group_numerator} must present in {.arg sample_info} (2nd column).")
  }

  if(!all(group_denominator %in% sample_info_groups)){
    cli::cli_abort("all values from {.arg group_denominator} must present in {.arg sample_info} (2nd column).")
  }


  # subset necessary columns from count.
  #The samples which are present in sample_information will only be processed further.

  count_data_select <- count_data %>%
    dplyr::select(!!column_geneid_quo, !!!sample_info_samples_quo)

  # display warning for the samples which won't be considered for analysis

  # // TO DO

  # in count_data convert NA to 0

  count_data_select <- .count_data_convert_NA_to_0(count_data = count_data_select)

  # remove rows having 0 in all samples discarded

  if(print_rows_all_zero){
    count_data_all_zero <- count_data_select %>%
      TidyWrappers::tbl_keep_rows_zero_all()

    if(nrow(count_data_all_zero) >0 ){
      rows_zero_all <- count_data_all_zero %>% dplyr::pull(1) %>% unique()

      cli::cli_text("Row{?s} {.emph {cli::col_red({rows_zero_all})}} have count 0 in all samples.
                {?That/Those} will be discarded.")
    }
  }

  # keep genes which have values non-zero in at least one sample.
  count_data_non_zero_in_one <- count_data_select %>%
    TidyWrappers::tbl_remove_rows_zero_all()

  count_data_filter <- .filter_rows_from_count_data(count_data_non_zero_in_one,
                                                    column_geneid = column_geneid,
                                                    sample_info = sample_info_data,
                                                    min_counts= min_counts,
                                                    min_replicates = min_replicates)

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

  # remove all combinations where numerator and denominator are same
  comb <- comb %>% dplyr::filter(numerator != denominator)

  if(nrow(comb) == 0){
    cli::cli_abort("{.arg group_numerator} and {.arg group_denominator} have at least a pair of different value.")
  }

  # summarize results
  cli::cli_alert_info("Summarizing DEG ...")

  # prepare a tibble where each row denotes to a combination of numerator and denominator.

  xx <- comb %>%

    ## add column comparisons
    dplyr::mutate(de_comparisons = stringr::str_c(.$numerator , .$denominator ,sep = "_VS_")) %>%

    ## add column norm counts
    dplyr::mutate(norm_counts = purrr::map2(.$numerator, .$denominator, ~ norm_counts[c(..1,..2)])) %>%


    # for each combination of numerator and denominator get deseq result. Results will be stored in a list column of tibble

    dplyr::mutate(dsr =  purrr::map2( .x = !!column_numerator,
                                      .y = !!column_denominator ,
                                      ~ DESeq2::results(dds, contrast = c("sample_groups" , ..1,..2)))) %>%

    # convert deseq result object into a tibble.

    dplyr::mutate(dsr_tibble = purrr::map(dsr , ~..1 %>% .dsr_to_tibble())) %>%

    # categorize DEG based on log2FC and pvalue

    dplyr::mutate(dsr_tibble_deg = purrr::map(dsr_tibble , ~ ..1 %>%
                                                .categorize_diff_genes(log2fc_cutoff = cutoff_lfc ,
                                                                       pval_cutoff = cutoff_pval,
                                                                       padj_cutoff = cutoff_padj ,
                                                                       regul_based_upon = regul_based_upon) )) %>%
    # summarize DEG

    dplyr::mutate(deg_summmary = purrr::map(dsr_tibble_deg , ~ ..1 %>%
                                              dplyr::group_by(regul) %>%
                                              dplyr::tally() ,.id = "cond"))


  ## use comparisons as names for list elements
  xx  %<>%  dplyr::mutate_if(is.list, ~rlang::set_names(. , xx$de_comparisons))

  # column rearrange

  xx %<>% dplyr::select(de_comparisons , tidyselect::everything())

  # assign names to each list in the tibble

  cli::cli_alert_success("Done.")

  ## set class

  class(xx) <- c("parcutils" , class(xx))
  return(xx)

}



#' Prepare a fold change matrix
#' @description This function returns a dataframe having first column gene names and subsequent columns are
#' fold change values for the comparisons passed through `sample_comparisons`.
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparisons a character vector denoting sample comparisons for which fold change values to be derived.
#' @param genes a character vector denoting gene names for which fold change values to be derived.
#'
#' @return a dataframe.
#' @export
#'
#' @examples
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t")
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#'
#' genes = parcutils::get_genes_by_regulation(x = res, sample_comparison = "treatment2_VS_control") %>% names()
#'
#' fc_df <- get_fold_change_matrix(x = res, sample_comparison = c("treatment2_VS_control" , "treatment1_VS_control"), genes = genes)
#'
#' print(fc_df)
#'
get_fold_change_matrix <- function(x , sample_comparisons , genes){

  # validate x

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))

  # validate sample_comparisons

  stopifnot("sample_comparisons must be a character vector." = is.character(sample_comparisons))

  # validate genes

  stopifnot("genes must be a character vector." = is.character(genes))

  # check if sample_comparisons present in x

  if(!any (x$de_comparisons %in% sample_comparisons)) {
    stop("None of values from 'sample_comparisons' present in the x.")
  }

  # get gene_id column name

  ## NOTE: There should be a way to get all user input parameters from the 'parcutils' object.

  gene_id_column <- x$dsr_tibble_deg[[1]] %>% colnames() %>% .[1]

  # subset data
  res <- x %>%
    dplyr::filter(de_comparisons %in% sample_comparisons ) %>%
    dplyr::pull(dsr_tibble_deg) %>%
    purrr::map(~ ..1 %>% dplyr::filter(!!rlang::sym(gene_id_column) %in% genes) %>% dplyr::select(1,3)) %>%
    dplyr::bind_rows(.id  = "comparisons") %>%
    tidyr::pivot_wider(id_cols = !!rlang::sym(gene_id_column) , names_from = comparisons, values_from = log2FoldChange)


  if(nrow(res) == 0){
    warning("None of the queried 'genes' are found in the 'x'.")
  }

  return(res)
}



#' Prepare a matrix of normalised gene expression values.
#' @description This function returns a dataframe having first column gene names and subsequent columns are
#' normalised gene expression values for the samples passed through argument `samples`. Internally it filters
#' the data from the  `norm_counts` of argument `x`. For more details on how `norm_counts` created refer to
#' the documentation of [parcutils::run_deseq_analysis()].
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting samples for which normalised gene expression values to be derived, Default NULL. If NULL it returns all samples in x
#' @param genes a character vector denoting gene names for which normalised gene expression values to be derived, Default NULL. If NULL it returns all samples in x.
#' @param summarise_replicates logical, default FALSE, indicating whether gene expression values summarised by mean or median between replicates.
#' @param summarise_method a character string either "mean" or "median" by which normalised gene expression values will be summarised between replicates.
#'
#' @return a dataframe.
#' @export
#'
#' @examples
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t")
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#'
#'
#' get_normalised_expression_matrix(x = res) %>% print()
#'
#' # summarise replicates by median
#'
#' get_normalised_expression_matrix(x = res ,summarise_replicates = T, summarise_method = "median") %>% print()
get_normalised_expression_matrix <- function(x , samples = NULL, genes = NULL, summarise_replicates = FALSE, summarise_method = "median" ){

  # validate x

  .validate_parcutils_obj(x)

  # validate genes

  stopifnot("genes must be a NULL or character vector." = is.character(genes) | is.null(genes))

  # validate summarise_replicates

  stopifnot("summarise_replicates must be a logical." = is.logical(summarise_replicates))

  # validate samples

  stopifnot("samples must be a NULL or character vector." = is.character(samples) | is.null(samples))

  # validate summarise_method

  match.arg(summarise_method , choices = c("mean" ,"median"))

  gene_id_column <- x$norm_counts[[1]][[1]] %>% colnames() %>% .[1]


  ## get all samples gene expression matrix.
  all_expr_mats <-  .get_all_named_expression_matrix(x)

  # Keep only user supplied samples
  if(is.null(samples)) {
    cols_selected_expr_mats <- all_expr_mats
  } else{
    cols_selected_expr_mats <- all_expr_mats[samples]
  }

  # make data long format.
  expr_mat_long <- cols_selected_expr_mats %>%
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
    if(!is.null(samples)){
      expr_mat_wide <- expr_mat_wide %>%
        dplyr::select(!!rlang::sym(column_gene_id),dplyr::all_of(samples))
    }

  } else{
    #  make it wide
    expr_mat_wide <- expr_mat_long %>%
      tidyr::pivot_wider(id_cols = !!rlang::sym(column_gene_id), values_from = "vals", names_from = "replicate" )

  }

  # filter by user supplied genes

  if(is.null(genes)) {
    ret <- expr_mat_wide
  }else{
    ret <- .filter_df_by_genes(df = expr_mat_wide, genes = genes)
  }

  return(ret)


}



#' Get genes based on their differential regulation
#' @details For a given differential comparison this function returns `up`, `down`, `both`, `other` and `all` genes.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparisons a character vector denoting  sample comparisons for which genes to be obtained.
#' @param regulation a character string, default \code{both}. Values can be one of the \code{up}, \code{down}, \code{both}, \code{other}, \code{all}.
#'  + `up` : returns all up regulated genes.
#'  + `down` : returns all down regulated genes.
#'  + `both` : returns all up and down regulated genes.
#'  + `other` : returns genes other than up and down regulated genes.
#'  + `all` : returns all genes.
#' @param simplify logical, default FALSE, if TRUE returns result in a dataframe format.
#' @return a list or dataframe.
#' @export
#'
#' @examples
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t")
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#'
#' # get both up and down regulated genes
#' get_genes_by_regulation(x = res, sample_comparisons = c("treatment1_VS_control")) %>% head()
#'
#' # get up genes only
#' get_genes_by_regulation(x = res, sample_comparisons = c("treatment1_VS_control") , regul = "up") %>% head()
#'
#' # get down genes only
#' get_genes_by_regulation(x = res, sample_comparisons = c("treatment1_VS_control") , regul = "down") %>% head()
#'
#' # get genes other than up and down
#' get_genes_by_regulation(x = res, sample_comparisons = c("treatment1_VS_control") , regul = "other") %>% head()
#'
#' # Simplify output for multiple sample comparisons
#' get_genes_by_regulation(x = res, sample_comparisons = res$de_comparisons, simplify = TRUE, regul= "up")
#'
#'
#' # get genesets by regulation. It uses sample comparison and regulation to name each output geneset.
#'
#' get_genesets_by_regulation(x = res, sample_comparisons = "treatment1_VS_control", regul = "both")
#'
get_genes_by_regulation <-  function(x, sample_comparisons , regulation = "both" , simplify = FALSE  ) {

  # validate x.
  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))

  # validate sample comparison.
  stopifnot("sample_comparisons must be a character vector" = is.character(sample_comparisons) & length(sample_comparisons) >= 1)

  # validate regulation.
  match.arg(regulation , choices = c("up","down" , "both", "other" ,"all"), several.ok = F)

  # validate simplify
  stopifnot("simplify must be a logical value" = is.logical(simplify))

  # if length of  sample_comparisons > 1, use this function recursively
  if(length(sample_comparisons) > 1){
    rslt <- purrr::map(sample_comparisons,
                       ~get_genes_by_regulation(x = x,sample_comparisons = ..1,regulation = regulation , simplify = simplify))
    names(rslt) <- sample_comparisons
    return(rslt)
  }

  # split by column regul.
  # NOTE: do not group by column name instead use index. The reason is because if the column name of column 'regul' change in future it will break this code.

  genes_by_comp <- x$dsr_tibble_deg[[sample_comparisons]] %>%
    dplyr::select(1, dplyr::last_col())

  genes_by_comp <- genes_by_comp %>%
    dplyr::group_by_at(2)

  grp_keys <- genes_by_comp %>%
    dplyr::group_keys() %>%
    dplyr::pull(1)

  genes_by_comp <- genes_by_comp %>%
    dplyr::group_split()
  names(genes_by_comp) <- grp_keys

  # Convert names and regulation to lower case before comparison
  names(genes_by_comp)  <- tolower(names(genes_by_comp))
  regulation = tolower(regulation)


  if(regulation == "both"){
    # select up and/or down
    rslt <- genes_by_comp[names(genes_by_comp) %in% c("up","down")]  %>% as.list() %>% dplyr::bind_rows()
  } else if(regulation == "all"){
    # select all
    rslt <- genes_by_comp  %>% as.list() %>% dplyr::bind_rows()
  } else{
    # select other
    rslt <- genes_by_comp[names(genes_by_comp) %in% regulation]  %>% as.list() %>% dplyr::bind_rows()
  }

  # rslt <- rslt[[1]] %>%
  #   purrr::set_names(rslt[[2]] %>% tolower())

  # convert factor
  if(nrow(rslt) == 0){
    rslt <- factor(levels = c("up","down","other"))
  } else{
    rslt <- factor(rslt[[2]] %>% tolower(), levels = c("up","down","other"))  %>%
      purrr::set_names(rslt[[1]])
  }

  if(simplify){
    rslt %<>% tibble::enframe() %>%
      dplyr::rename(c("regul" = "value" , "genes" = "name")) %>%
      dplyr::mutate(sample_comparisons = sample_comparisons) %>%
      dplyr::select(dplyr::all_of(c("sample_comparisons", "genes", "regul")))
  }
  return(rslt)

}

#' @rdname get_genes_by_regulation
#' @export
get_genesets_by_regulation <- function(x, sample_comparisons, regulation = "both" ){


  genes <- get_genes_by_regulation(x, sample_comparisons = sample_comparisons,
                                   regulation = regulation,
                                   simplify = TRUE)

  gene_sets <- dplyr::bind_rows(genes) %>%
    dplyr::mutate(gene_set_name = stringr::str_c(sample_comparisons, regul , sep = "_")) %>%
    named_group_split(gene_set_name) %>%
    purrr::map(~ ..1 %>% dplyr::pull(genes))

  return(gene_sets)
}






#' Generate upset plots for differently expressed genes between comparisons.
#' @description  For a given set of DEG comparisons, the functions returns [UpSetR::upset()] plots for up and down genes between any 2 comparisons.
#' For each upset plot generated function also returns interaction between gene sets in form of dataframe.
#' @param x x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparisons a character vector denoting  sample comparisons between upset plot to be generated.
#' @param color_up a character string denoting a valid color code for bars in upset plot for up regulated genes.
#' @param color_down a character string denoting a valid color code for bars in upset plot for down regulated genes.
#'
#' @return an object of named list where each element is a list of two - 1)  an upset plots  [UpSetR::upset()] and their intersections in form of tibble.
#' @export
#' @examples
#'
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t")
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#'
#' yy <- plot_deg_upsets(x= res, sample_comparisons = c("treatment1_VS_control","treatment2_VS_control"))
#'
#' yy[[1]]$upset_plot %>% print()
#'
#' yy[[1]]$upset_intersects %>% print()
#'
plot_deg_upsets <- function(x, sample_comparisons, color_up = "#b30000", color_down = "#006d2c"){

  # x <- dd
  # sample_comparisons = c("shRNA1_VS_shCTRL", "INFy_shRNA1_VS_INFy_shCTRL" , "INFy_shCTRL_VS_shCTRL")

  # validate x

  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))


  # validate sample_comparisons

  stopifnot("sample_comparisons must be a character vector of length 2." = is.character(sample_comparisons) & length(sample_comparisons) > 1)

  # all sample comparisons must present in x

  if(!all(sample_comparisons %in% x$de_comparisons) ){
    not_present <- sample_comparisons[!(sample_comparisons %in% x$de_comparisons)]
    stop(glue::glue("Value {not_present} is not present in the x. Check x$de_comparisons to access all choices."))
  }

  # validate colors

  stopifnot("color_up and color_down must be a character string denoting a valid color code." = is.character(color_up) & is.character(color_up) & length(color_up)==1, length(color_down)==1)

  # generate all combinations of 2 from sample_comparisons.
  sample_comparisons_all_comb <- combn(sample_comparisons, 2)  %>%
    as.data.frame() %>%
    purrr::map( ~..1)


  ## set names
  sample_comparisons_all_comb <- purrr::set_names(sample_comparisons_all_comb , ~purrr::map_chr(sample_comparisons_all_comb , ~paste0(..1, collapse = "_AND_")))


  all_upsets <- purrr::map(sample_comparisons_all_comb, ~ .piarwise_upset(x = x,
                                                                          sample_comparison = ..1,
                                                                          color_up = color_up ,
                                                                          color_down = color_down) )

  return(all_upsets)
}



#' Generate a heatmap of normalised gene expression values, z-score or log2Fold-change values.
#' @description Heatmap is a common tool to show gene expression pattern across samples in RNA-seq experiments.
#' During data exploration, it is a common practice to generate several heatmaps to identify interesting
#' patterns of gene expression across samples. Most heatmap generating tools require data in a tabular format.
#' However, to prepare such a data requires significant data wrangling such as sub-setting
#' genes (rows) and samples (columns) from multiple R objects. Such data wrangling creates redundant, and less readable code
#' which increases the chances of error. Furthermore, because of less readability using the same code for another analysis
#' requires code cleaning and code rewriting making code less reusable.
#' Situations complicates even more when RNA-seq studies consist of several samples and several comparisons.
#' This function cut-downs several steps of data wrangling to create a heatmap of gene expression, z-score or
#' log2 fold-change. Given an object of the class \code{parcutils}, sample names or sample comparisons, genes to show in the heatmap
#' and several other arguments, the heatmap can be created quickly. The output heatmap is an output of the function [ComplexHeatmap::Heatmap()]
#' which then can be used with other functions of the [ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/) package.
#'
#' @param x an abject of the class \code{parcutils}. This is can be created using the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting sample names to use in the heatmap.
#' @param genes a character vector denoting genes to use in the heatmap.
#' @param repair_genes  logical, default \code{FALSE}, indicating whether to repair gene names. See details.
#' @param convert_log2 logical, default \code{FALSE}, indicating whether to log2 transform gene expression values.
#' @param color_default logical, default \code{TRUE}, indicating whether to use default heatmap colors.
#' @param col an output of [circlize::colorRamp2()], default NULL.
#' @param convert_zscore logical, default \code{TRUE}, indicating whether to convert gene expression values in to the z-score or not. see details.
#' @param summarise_replicates logical, default \code{TRUE}, indicating whether to summarise values for each gene across replicates.
#' @param summarise_method a character string, either mean or median.
#' @param show_row_names logical, default \code{FALSE}, indicating whether to show row names in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param cluster_rows logical, default \code{TRUE}, indicating whether to cluster rows in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param show_row_dend logical, default \code{TRUE}, indicating whether to show dendrogram in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param row_names_font_size a numeric value, default 10, indicating size of row names in the heatmap.
#' @param show_column_names logical, default \code{TRUE}, indicating whether to show column names in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param cluster_columns logical, default \code{TRUE}, indicating whether to cluster columns in the heatmap or not.
#' Internally this argument is passed to the function [ComplexHeatmap::Heatmap()].
#' @param show_heatmap_legend logical, default \code{TRUE}, indicating whether to show heatmap legend or not.
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
#'
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t")
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#'
#' genes = parcutils::get_genes_by_regulation(x = res, sample_comparison = "treatment2_VS_control" , "both") %>% names()
#' get_gene_expression_heatmap(x = res, samples = c("control" ,"treatment1" , "treatment2"), genes = genes)
#'
#' # plot raw expression values
#'
#' get_gene_expression_heatmap(x = res, samples = c("control" ,"treatment1" , "treatment2"),
#' genes = genes, convert_zscore = FALSE)
#'
#' # plot log2 expression values
#'
#' get_gene_expression_heatmap(x = res, samples = c("control" ,"treatment1" , "treatment2"),
#' genes = genes, convert_zscore = FALSE,convert_log2 = TRUE)
#'
#' # plot all replicates
#'
#' get_gene_expression_heatmap(x = res, samples = c("control" ,"treatment1" , "treatment2"),
#' genes = genes, convert_zscore = TRUE, summarise_replicates = FALSE)
#'
#' # show gene names
#'
#' get_gene_expression_heatmap(x = res, samples = c("control" ,"treatment1" , "treatment2"),
#' genes = genes[1:10], convert_zscore = TRUE, summarise_replicates = FALSE , show_row_names = TRUE)
#'
#' # repair gene names
#'
#'get_gene_expression_heatmap(x = res, samples = c("control" ,"treatment1" , "treatment2"),
#' genes = genes[1:10], convert_zscore = TRUE, summarise_replicates = FALSE, show_row_names = TRUE, repair_genes = TRUE)
#'
#"
#' # fold change heatmap
#'
#' get_fold_change_heatmap(x = res , sample_comparisons = c("treatment2_VS_control" ,"treatment1_VS_control") ,
#' genes = genes , cluster_columns = FALSE , name = "Log2FC")
#'
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
  all_expr_mat <-  .get_all_named_expression_matrix(x)

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

  expr_mat_wide <- .filter_df_by_genes(df = expr_mat_wide, genes = genes)

  # remove row if all are NA.
  value_na <- expr_mat_wide  %>%  TidyWrappers::tbl_keep_rows_NA_any() %>%
    dplyr::pull(!!rlang::sym(column_gene_id))
  if(length(value_na) > 0){
    value_na <- paste0(value_na, collapse = ",")
    cli::cli_alert_warning("Gene{?s} with value NA - {.emph {cli::col_red({value_na})}} - {?is/are} removed from heatmap.")
    expr_mat_wide <- expr_mat_wide %>% TidyWrappers::tbl_remove_rows_NA_any()
  }

  #fix colors
  if(color_default){
    col = .fix_hm_colors(expr_mat_wide)
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
    col = .fix_hm_colors(fc_df)
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







#' Get data from a heatmap in the same order.
#' @description extract genes or gene clusters from the heatmap in same order.
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


#' Compare log2 fold change between two sample comparisons.
#' @description This function generates a scatter plot of log2 fold change values for two different comparisons.
#'
#' @param x an object of class parcutils.
#' @param sample_comparisons a character vector of length 2 denoting sample comparisons to plot.
#' @param labels a character vector of genes to label. Default NULL, show all genes.
#' @param point_size a numeric, default 2, denoting size of the points.
#' @param label_size a numeric, default 2, denoting size of the labels.
#' @param col_up a character string, default `#a40000`, a valid color code for common up regulated genes.
#' @param col_down a character string, default `#16317d`, a valid color code for common down regulated genes.
#' @param color_label a character string one of the "both", #both_down, or "both_up". Default `both`.
#' @param repair_genes a logical, default `TRUE`, denotes whether to repair gene names or not,
#' If `TRUE` string prior to `:` will be removed from the gene names.
#'
#' @return an object of ggplot2.
#' @export
#' @examples
#' count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
#' count_data <- readr::read_delim(count_file, delim = "\t")
#'
#'sample_info <- count_data %>% colnames() %>% .[-1]  %>%
#'  tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )
#'
#'
#'res <- run_deseq_analysis(counts = count_data ,
#'                          sample_info = sample_info,
#'                          column_geneid = "gene_id" ,
#'                          group_numerator = c("treatment1", "treatment2") ,
#'                          group_denominator = c("control"))
#' # show common up and common down
#' get_fold_change_scatter_plot(x = res,
#' sample_comparisons = c("treatment1_VS_control",
#' "treatment2_VS_control"),label_size = 3)
#'
#' # show common up
#'
#' get_fold_change_scatter_plot(x = res,
#' sample_comparisons = c("treatment1_VS_control",
#' "treatment2_VS_control"),
#' color_label = "both_up",label_size = 4)
#'
#'  # show common down
#' get_fold_change_scatter_plot(x = res,
#' sample_comparisons = c("treatment1_VS_control",
#' "treatment2_VS_control"),
#' color_label = "both_down",label_size = 4, point_size = 4)
#'
get_fold_change_scatter_plot <- function(x,
                                         sample_comparisons,
                                         labels = NULL,
                                         point_size = 2,
                                         label_size = 2,
                                         color_label = "both", #both_down, "both_up"
                                         col_up = "#a40000",
                                         col_down = "#16317d",
                                         repair_genes = T
){

  # validate arguments

  .validate_parcutils_obj(x)

  # sample_comparisons
  if(!(is.character(sample_comparisons) & length(sample_comparisons) == 2)){
    cli::cli_abort("{.arg sample_comparisons} must be a {.cls character} vector of length 2.")
  }

  #sample_comparisons must present in x$de_comparisons

  .is_sample_comparsions_present_in_obj_parcutils(x, sample_comparisons = sample_comparisons)

  # labels
  if(!(is.null(labels) | is.character(labels))){
    cli::cli_abort("{.arg labels} must be a {.cls NULL} or a {.cls character} vector.")
  }

  # point_size
  if(!(is.numeric(point_size) & length(point_size) == 1)){
    cli::cli_abort("{.arg point_size} must be a {.cls numeric} value of length 1.")
  }

  # label_size
  if(!(is.numeric(label_size) & length(label_size) == 1)){
    cli::cli_abort("{.arg label_size} must be a {.cls numeric} value of length 1.")
  }
  # color_label
  match.arg(color_label , choices = c("both","both_up" ,"both_down"))

  #col_up
  if(!(is.character(col_up) & length(col_up) == 1)){
    cli::cli_abort("{.arg col_up} must be a {.cls character} value denoting a valid color name.")
  }
  #col_down
  if(!(is.character(col_down) & length(col_down) == 1)){
    cli::cli_abort("{.arg col_down} must be a {.cls character} value denoting a valid color name.")
  }
  #repair_genes
  if(!is.logical(repair_genes)){
    cli::cli_abort("{.arg repair_genes} must be a {.cls logical} value.")
  }

  # check presence of labels in x
  .is_genes_present_in_obj_parcutils(x = x, genes = labels)

  # prepare gene groups
  gene_groups <- .get_gene_groups_for_fc_scatter_plot(x,
                                                      sample_comparisons = sample_comparisons,
                                                      color_label = color_label
  )

  gene_groups_col_names <- gene_groups %>% colnames()

  # prepare color vec
  if(color_label == "both_up"){
    color_group <- c("both_up" = col_up,
                     "other" = "grey")

  }else if(color_label == "both_down"){
    color_group <- c("both_down" = col_down,
                     "other" = "grey")

  }else{
    color_group <- c("both_up" = col_up,
                     "both_down" = col_down,
                     "other" = "grey")

  }

  # get all genes from x to plot
  genes_for_plot <- .get_all_expressed_genes(x = x) %>% names()

  # get fold change values for the required comparisons.
  fc_data <- get_fold_change_matrix(x = x,
                                    sample_comparisons = sample_comparisons,
                                    genes = genes_for_plot)

  fc_data_col_names <- fc_data %>% colnames()

  # join fc values with gene groups
  for_plot <- fc_data %>%
    dplyr::left_join(gene_groups,
                     by = purrr::set_names(gene_groups_col_names[1],fc_data_col_names[1]))


  # plot
  gp <- for_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = !!rlang::sym(sample_comparisons[1]),
                                 y = !!rlang::sym(sample_comparisons[2]))) +
    ggplot2::geom_point(aes(col = group)) +

    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = color_group)

  # prepare data for labels
  if(!is.null(labels)){
    data_for_labels  <- for_plot %>%
      dplyr::filter(!!rlang::sym(fc_data_col_names[1]) %in% labels)
  } else{
    data_for_labels <- for_plot %>%
      dplyr::filter(group != "other")
  }

  # repair genes
  if(repair_genes){
    data_for_labels %<>%
      dplyr::mutate(gene_id = stringr::str_replace(string = !!rlang::sym(fc_data_col_names[1]),
                                                   pattern = ".*:",
                                                   replacement = ""
      ))
  }

  gp <- gp + ggrepel::geom_text_repel(data =  data_for_labels,
                                      aes(label = !!rlang::sym(fc_data_col_names[1])),
                                      size = label_size)
  # add suffix 'Log2FC' to axis labels
  suffix <- "Log2FC"
  gp$labels$x <- glue::glue("{suffix}\n({gp$labels$x})")
  gp$labels$y <- glue::glue("{suffix}\n({gp$labels$y})")
  return(gp)

}







