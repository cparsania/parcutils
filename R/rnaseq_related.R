#' @title Perform differential expression analysis using `DESeq2`
#' @description This is a wrapper function build upon `DESeq2::DESeq()` and `DESeq2::DESeqResults()`
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

  ## prepare annotation ----
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


#' Convert DESeq result object in to a tibble.
#'
#' @param x \link{DESeqResults} object
#'
#' @return a data frame
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate_if
#' @importFrom tibble as_tibble
#' @examples
#' \dontrun{
#'
#' }
#'
dsr_to_tibble <- function(x){

  x %>% as.data.frame() %>%
    tibble::rownames_to_column("GeneID") %>%
    tibble::as_tibble() #%>%
    #dplyr::mutate_if(is.numeric, ~ round(..1, 3)) ## round by 3 digits
}





##@@@@@@@@
## filter gff ----
## filter gff
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
#'
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
#'
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



### export and visualizations ----

# # map annotations
#
# dsr_tibble <- purrr::map(dsr_tibble, ~ .x %>% dplyr::left_join(gtf_annot_gene , c("Geneid" = "gene_id")))
#
# # summarise diff genes
#
# diff_genes_count <- dsr_tibble %>%
#   purrr::map_df(~ ..1 %>% dplyr::group_by(regul) %>%
#                   dplyr::tally() ,.id = "cond")
#
# # export diff data to file
#
# diff_data_paths <- stringr::str_c("./","deseq_out_",outfile_prefix ,"_", names(dsr_tibble) ,".txt")
# purrr::pwalk(list(dsr_tibble, diff_data_paths) , readr:::write_delim , delim = "\t" )
#
#
# #purrr::map2(comp, dsr_tibble , ~ readr:::write_delim(..2 , file = glue::glue("./diff_data_{..1}.txt"),delim = "\t"))
#
#
#
# ## export normalised counts ---
#
# normalized_counts <-
#   DESeq2::counts(dds , normalized = T) %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("Geneid") %>%
#   tibble::as_tibble()
#
# normalized_counts_log2 <- normalized_counts %>%
#   TidyWrappers::tbl_convert_log2(frac = 1)
#
# norm_count_file_path <- stringr::str_c("./",outfile_prefix,"_","normalised_count_log2.txt")
#
# readr::write_delim(normalized_counts_log2 , norm_count_file_path ,delim = "\t")
#

#
# # diff genes bar plot
#
# bar_plot_data <-
#   purrr::map_df(dsr_tibble, ~ ..1 %>% dplyr::filter(regul != "other") %>%
#                   dplyr::select(gene_name, regul) %>%
#                   dplyr::group_by(regul) %>%
#                   dplyr::tally() ,.id = "cond")
#
# bar_plot <- bar_plot_data %>%
#   ggplot2::ggplot() +
#   geom_bar(aes(x = regul, y = n, fill = regul), stat = "identity") +
#   ggplot2::facet_wrap(~cond) +
#   theme_bw() +
#   theme(text = element_text(size = 15)) +
#   scale_fill_manual(values = c("Down" = "green3" , "Up" ="red2")) +
#   ylab("No. of Genes") +
#   xlab("Type of regulation")
#
# bar_plot_paths <- stringr::str_c(outfile_prefix,"_","diff_count_bar_plot", ".pdf")
# ggsave(filename = bar_plot_paths, bar_plot, height = 5, width = 7 , units = "in")
#
#
# # normalised count long format data
#
# normalized_counts_long <- normalized_counts_log2 %>%
#   tidyr::gather(key = samples, value = value , -Geneid)
#
# # add sample information to long data
#
# normalized_counts_long %<>% dplyr::left_join(sampleinfo, c("samples" = "id"))
#
# # PCA plot
#
# pca_input <- normalized_counts_log2 %>%
#   tidyr::gather(key = sample ,value = value , -Geneid) %>%
#   tidyr::spread(Geneid, value) %>% as.data.frame() %>%
#   tibble::column_to_rownames("sample")
#
# pr_comps <- prcomp(pca_input,scale = T)
#
# pr_comps_tbl <- pr_comps %>%
#   .$x %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column("samples") %>%
#   tibble::as_tibble()
#
# # get prop of variance
#
# pc_prop_of_var <- summary(pr_comps)$importance %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "feature_type") %>%
#   tidyr::gather(pc , value , -feature_type) %>%
#   dplyr::filter(feature_type == "Proportion of Variance") %>%
#   dplyr::select(-feature_type) %>%
#   dplyr::mutate(value = round(value * 100 ,1)) %>% ## convert to percentage
#   dplyr::mutate(with_var = paste(pc ," (" , value ,"%",")" , sep = "")) %>%
#   dplyr::select(-value)
#
# # return named vector
#
# pc_prop_of_var <- pc_prop_of_var %>%
#   dplyr::pull(with_var) %>%
#   rlang::set_names( pc_prop_of_var$pc)
#
# # add sample info
#
# pr_comps_tbl %<>% dplyr::left_join(sampleinfo , by = c("samples" = "id"))
#
# # select which PCs to plot
#
# x_pc = quo(PC1)
# y_pc = quo(PC2)
#
# # pca plot
#
# pca_plot <- pr_comps_tbl %>%
#   ggplot2::ggplot(aes(x = !!x_pc, y = !!y_pc)) +
#   geom_point(aes(col = cond),size = 5) +
#   theme_bw()+
#   theme(text = element_text(size = 15,colour = "black"))+
#   xlab(pc_prop_of_var[[rlang::as_label(x_pc)]]) +
#   ylab(pc_prop_of_var[[rlang::as_label(y_pc)]])
#
# # save pca plot
#
# pca_plot_paths <- stringr::str_c(outfile_prefix,"_","pca_plot", ".pdf")
# purrr::pwalk(list(pca_plot_paths,list(pca_plot)), ggsave, width = 6, height = 5, units = "in")
#
# ggsave(filename = glue::glue("{outfile_prefix}_PCAplot.pdf") , width = 6 , height = 5, units = "in")
#
# # samplewise box plot
#
# box_plots <- normalized_counts_long %>%
#   ggplot2::ggplot() +
#   geom_boxplot(aes(x = samples, y = value, fill =  cond)) +
#   theme_bw()+
#   theme(text = element_text(size = 15) , axis.text.x  = element_text(angle = 90)) +
#   ylab("Log2(normalised_count + 1)")
#
# box_plot_paths <- stringr::str_c(outfile_prefix,"_","box_plot", ".pdf")
#
# purrr::pwalk(list(box_plot_paths,list(box_plots)), ggsave, width = 6.5, height = 6, units = "in")
#
# # sample wise density plot
#
# density_plots <- normalized_counts_long %>%
#   ggplot2::ggplot() +
#   geom_density(aes(x = value, fill =  cond)) +
#   facet_wrap(~samples)+
#   theme_bw() +
#   theme(text = element_text(size = 15)) +
#   xlab("Log2(normalised_count + 1)")
#
# # save density plot
#
# density_plot_paths <- stringr::str_c(outfile_prefix,"_","density_plot", ".pdf")
#
# purrr::pwalk(list(density_plot_paths,list(density_plots)), ggsave, width = 6.5, height = 6, units = "in")
#
#
# # pairwise corr
#
# coldata_grpd <-
#   coldata %>%
#   tibble::rownames_to_column("id") %>%
#   dplyr::group_by(cond) %>%
#   tidyr::nest() %>%
#   dplyr::mutate(data = purrr::map(data, ~ ..1 %>% dplyr::pull(1)))
#
# coldata_lst <-
#   coldata_grpd %>%
#   dplyr::pull(2)
#
# names(coldata_lst) <- coldata_grpd %>% dplyr::pull(1)
#
# replicate_corr_plts <- purrr::map(coldata_lst , ~ GGally::ggpairs(normalized_counts_log2,
#                                                                   columns = ..1) + theme_bw() +
#                                     theme(text = element_text(size = 15)))
#
# # save pairwise corr plots
#
# corr_plot_paths <- stringr::str_c(outfile_prefix,"_",names(replicate_corr_plts),"_corr_plot", ".pdf")
#
# purrr::pwalk(list(corr_plot_paths, replicate_corr_plts), ggsave, width = 6, height = 6, units = "in")
#
# # volcano plot
#
# # plot volcano for all condition
#
# log2fc_cutoff = 1.3
# pval_cutoff = 0.05
#
# volcano_x_intercept = log2fc_cutoff  # log2fc
# volcano_y_intercept = pval_cutoff  # pvalue
#
#
# # prepare volcano plots
#
# volcano_plts <- purrr::map(names(dsr_tibble) , ~ dsr_tibble[[..1]] %>%
#                              ggplot(aes(x = log2FoldChange , y = -log10(pvalue))) +
#                              geom_point(aes(col = type),size =2, alpha = 0.5) +
#                              geom_vline(xintercept = c(-log2fc_cutoff,log2fc_cutoff)) +
#                              geom_hline(yintercept = c(-log10(pval_cutoff))) +
#                              # add lables
#                              ggrepel::geom_text_repel(data = dsr_tibble[[..1]] %>%
#                                                         dplyr::filter((log2FoldChange <= -1.8 | log2FoldChange >= 1.8)
#                                                                       & dplyr::between(-log10(pvalue) , 1.5,3)),
#                                                       aes(label = gene_name), max.overlaps = 20,size = 3) +
#
#                              theme_bw() +
#                              theme(text = element_text(size = 15)) +
#                              ggtitle(..1))
#
#
# names(volcano_plts) <- names(dsr_tibble)
#
#
# # save volcano plots
#
# volcano_paths  <- stringr::str_c(outfile_prefix,"_",names(volcano_plts),"_volcano_plot", ".pdf")
#
# purrr::pwalk(list(volcano_paths,volcano_plts) , ggsave, width = 8,height = 7, units = "in")
#
# # MA plot
#
# ma_fc_threshold <- 3
# ma_padj_threshold <- 0.01
#
# # prepare data for MA plot
#
# for_ma <- purrr::map(dsr_tibble,
#                      ~ ..1 %>% dplyr::mutate(log2FoldChange = dplyr::case_when(log2FoldChange >= ma_fc_threshold ~ ma_fc_threshold ,
#                                                                                log2FoldChange <= -1 * (ma_fc_threshold) ~ -1 * (ma_fc_threshold),
#                                                                                TRUE ~ log2FoldChange)) %>%
#                        dplyr::mutate(is_sig_by_padj = dplyr::case_when(padj <= ma_padj_threshold ~ TRUE,TRUE ~FALSE))
#
# )
#
# # generate MA plot
#
# ma_plots <- purrr::map(names(for_ma) , ~ for_ma[[..1]] %>%
#                          ggplot(aes(x = baseMean , y = log2FoldChange)) +
#                          geom_point(size = 0.5, aes(col = is_sig_by_padj)) +
#                          ggtitle(..1) +
#                          scale_x_log10() +
#                          geom_hline(yintercept = 0 , size =1 , col = "grey40") +
#                          theme_bw() +
#                          scale_color_manual(values =  c("TRUE" ="blue" , "FALSE" = "grey")) +
#                          theme(text = element_text(size = 15))+
#                          labs(color = "padj <= 0.01")
# )
#
# names(ma_plots) <- names(for_ma)
#
# # save MA plots
#
# ma_plot_paths <- stringr::str_c(outfile_prefix,"_",names(ma_plots),"_MA_plot", ".pdf")
#
# purrr::pwalk(list(ma_plot_paths, ma_plots), ggsave, width = 7,height = 6, path = ".")
#
#
# # upset plot for diff genes
#
# # upset plot requires a list containing character string
#
# # prepare a list of up and down genes from each diff comparisons
#
# diff_genes <- purrr::map_df(dsr_tibble , ~ ..1 %>%
#                               dplyr::filter(regul %in% c("Up","Down")) %>%
#                               dplyr::select(Geneid,gene_name, regul, log2FoldChange, pvalue, padj) ,
#                             .id = "cond")
#
# diff_genes_grpd <- diff_genes %>%
#   dplyr::mutate(group_var  = stringr::str_c(cond, "_",regul)) %>%
#   dplyr::select(gene_name, group_var) %>%
#   dplyr::group_by(group_var)
#
# grps_keys <- dplyr::group_keys(diff_genes_grpd) %>% dplyr::pull(1)
#
# upset_data <- diff_genes_grpd %>%
#   dplyr::group_split() %>%
#   purrr::map(~ ..1 %>% dplyr::pull(1))
#
# names(upset_data) <- grps_keys
#
# us_plot <- UpSetR::upset(UpSetR::fromList(upset_data) ,
#                          order.by =  "degree",
#                          group.by = "degree",
#                          nsets = length(upset_data),
#                          sets = names(upset_data) %>% rev(),
#                          keep.order=T,
#                          text.scale=1.5,
#                          sets.x.label = "No. of genes")
#
# # save upset plot
# us_plot_paths <- stringr::str_c(outfile_prefix,"_","upset_plot", ".pdf")
# pdf(file = us_plot_paths,width = 6 , height = 5)
# us_plot
# dev.off()
#
# # GO analysis
#
# go_input <- diff_genes %>%
#   dplyr::mutate(group_var = stringr::str_c(cond , "_",regul)) %>%
#   dplyr::select(group_var , Geneid) %>%
#   dplyr::group_by(group_var)
#
#
# group_keys <- dplyr::group_keys(go_input)
# go_input_lst <- go_input %>%
#   dplyr::group_split() %>%
#   purrr::map(~ ..1 %>% dplyr::pull(Geneid))
#
# names(go_input_lst) <- group_keys$group_var
#
# lengths(go_input_lst)
#
# ego1 <- clusterProfiler::enrichGO(gene  = go_input_lst$cond_treatment_vs_control_Down,
#                                   OrgDb         = org.Hs.eg.db,
#                                   keyType       = 'ENSEMBL',
#                                   ont           = "BP",
#                                   pAdjustMethod = "BH",
#                                   pvalueCutoff  = 0.01,
#                                   qvalueCutoff  = 0.05)
#
# p1  <- clusterProfiler::emapplot(enrichplot::pairwise_termsim(ego1) ,
#                                  showCategory = 30,
#                                  color = "p.adjust",
#                                  layout = "kk") + ggtitle("control_vs_treatmet_down")
#
# ggsave(p1, filename = "control_vs_treatment_Down.pdf" , height = 8, width = 8 , units = "in")
# ego2 <- enrichGO(gene         = go_input_lst$cond_treatment_vs_control_Up,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENSEMBL',
#                  ont           = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)
#
# p2  <- emapplot(pairwise_termsim(ego2),
#                 showCategory = 20,
#                 color = "p.adjust",
#                 layout = "kk") + ggtitle("control_vs_treatment_Up")
#
# ggsave(p2, filename = "cond_5KD_vs_control_Up.pdf" , height = 8, width = 8 , units = "in")
#
# ego3 <- enrichGO(gene         = go_input_lst$cond_6KD_vs_control_Down,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENSEMBL',
#                  ont           = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)
#
# p3  <- emapplot(pairwise_termsim(ego3),
#                 showCategory = 30,
#                 color = "p.adjust",
#                 layout = "kk") + ggtitle("cond_6KD_vs_control_Down")
#
# ggsave(p3, filename = "cond_6KD_vs_control_Down.pdf" , height = 8, width = 8 , units = "in")
#
# ego4 <- enrichGO(gene         = go_input_lst$cond_6KD_vs_control_Up,
#                  OrgDb         = org.Hs.eg.db,
#                  keyType       = 'ENSEMBL',
#                  ont           = "BP",
#                  pAdjustMethod = "BH",
#                  pvalueCutoff  = 0.01,
#                  qvalueCutoff  = 0.05)
#
# p4  <- emapplot(pairwise_termsim(ego4),
#                 showCategory = 20,
#                 color = "p.adjust",
#                 layout = "kk") + ggtitle("cond_6KD_vs_control_Up")
#
# ggsave(p4, filename = "cond_6KD_vs_control_Up.pdf" , height = 8, width = 8 , units = "in")
#
#



