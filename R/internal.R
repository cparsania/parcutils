
#' Filter expression matrix (dataframe) by genes (first column).
#'
#' @param df a gene expression data frame
#' @param genes a character vector of genes to be filtered from df.
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
.filter_df_by_genes <- function(df, genes){

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


#' Get named list of gene expression matrix for all samples in \code{x}.
#'
#' @param x an object of class parcutils.
#'
#' @return named list of gene expression matrix.
#' @export
#' @keywords internal
.get_all_named_expression_matrix <- function(x){

  .validate_parcutils_obj(x)

  all_sample_names <- x$norm_counts %>%
    purrr::flatten() %>% names()

  unique_index <- purrr::map_dbl(all_sample_names %>% unique() , ~ which(..1 == all_sample_names)[1])

  ret <- x$norm_counts %>% purrr::flatten() %>% magrittr::extract(unique_index)
  return(ret)

}


#' Validate if an object is of class 'parcutils'.
#'
#' @param x an object
#'
#' @return error
#' @export
#' @keywords internal
.validate_parcutils_obj <- function(x){
  stopifnot("x must be an object of class 'parcutils'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils"))
}


#' Given a matrix generate a color scale for a heatmap.
#'
#' @param hm_matrix
#'
#' @return color for [ComplexHeatmap::Heatmap()].
#' @export
#' @keywords internal
.fix_hm_colors <- function(hm_matrix){

  stopifnot(is.data.frame(hm_matrix))

  col_choices <- c("min" = MetBrewer::met.brewer(name = "Austria")[2],
                   "max" = MetBrewer::met.brewer(name = "Austria")[1])

  max_val <-  ceiling(hm_matrix %>% dplyr::select_if(is.numeric) %>% max())
  min_val <-  floor(hm_matrix %>% dplyr::select_if(is.numeric) %>% min() )
  mid_val <- (max_val + min_val) / 2

  col = circlize::colorRamp2(c(min_val , mid_val, max_val), c(col_choices["min"], "white", col_choices["max"]))


}

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
#'
.piarwise_upset <- function(x, sample_comparison , color_up = "#b30000", color_down = "#006d2c",... ){

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



#' @title  Categorize diff genes in to `up` and `down` based on log2FC, p-value and p-adj values.
#' @description Usually fold change 1.5 (log2FC ~0.6) and statistical significance defined by p and/or p.adj values
#' are considered to mark gene as differential expressed. However, these cutoffs change based on data quality.
#' This function allows changing these cutoffs and mark genes differently expressed.
#' While using this package, most of the time user do not need to call this function explicitly as it is internally called by
#' [parcutils::run_deseq_analysis()].
#'
#' @param dsr_tibble  a dataframe obtained from DESeqResult object or an object of DESeqResult.
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
#' @keywords internal
.categorize_diff_genes <-  function(dsr_tibble,
                                   log2fc_cutoff =  1,
                                   pval_cutoff = 0.05,
                                   padj_cutoff = 0.01,
                                   add_column_regul = TRUE,
                                   regul_based_upon = 1
){

  # convert DESeqResult object to tibble

  if(any(is(dsr_tibble) %in% "DESeqResults")) {
    dsr_tibble %<>% .dsr_to_tibble()
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



#' Filter rows from a GFF file.
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
.filter_gff <- function(gtf_file ,
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


#' Convert DESeq result object in to a tibble.
#'
#' @param x DEseq2 result object. Can be generated by [DESeq2::DESeqResults()].
#'
#' @return a data frame.
#' @export
#' @keywords internal
#'
.dsr_to_tibble <- function(x, .col_gene_id = "gene_id"){

  x %>% as.data.frame() %>%
    tibble::rownames_to_column(.col_gene_id) %>%
    tibble::as_tibble() #%>%
  #dplyr::mutate_if(is.numeric, ~ round(..1, 3)) ## round by 3 digits
}



#' Group replicates by samples.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#'
#' @return a tbl.
#' @export
#' @keywords internal
.group_replicates_by_sample <- function(x){

  .validate_parcutils_obj(x)

  sample_info <- purrr::map(.get_all_named_expression_matrix(x), ~ ..1 %>% colnames() %>% .[-1]) %>%
    tibble::tibble(groups = names(.) , samples = .) %>%
    tidyr::unnest(cols = "samples")

  return(sample_info)
}



#' Get list of all expressed genes from \code{x}
#'
#' @param x an object of class 'parcutils'
#'
#' @return a character vector.
#' @export
#' @keywords internal
.get_all_expressed_genes <- function(x){
  .validate_parcutils_obj(x)

  all_genes <- parcutils::get_genes_by_regulation(x = x,
                                                  regulation = "all",
                                                  sample_comparisons = x$comp[1])

  return(all_genes)
}



#' Remove items from the list if no enriched terms found.
#'
#' @param x a list containing output of GO enrichment. Each element of the list is an output of  [clusterProfiler::enrichGO()].
#'
#' @return a list of containing output of GO enrichment.
#' @export
#' @keywords internal
.keep_only_enriched_go <- function(x){

  safely_enrich <- purrr::safely(enrichplot::pairwise_termsim , otherwise = NULL)

  ## remove those which has no enriched term

  go_enrichment_result_only_enriched <- x %>%
    purrr::map( ~ safely_enrich(..1)) %>%
    purrr::transpose() %>%
    purrr::simplify_all()

  go_enrichment_result_only_enriched <-
    go_enrichment_result_only_enriched$result %>%
    purrr::discard(is.null)

}


#' Remove highly similar GO.
#'
#' @param x a list containing output of GO enrichment. Each element of the list is an output of  [clusterProfiler::enrichGO()].
#' @param similarity_cutoff a value between 0 and 1. Default 0.9.
#'
#' @return  a list containing output of GO enrichment.
#' @export
#' @keywords internal
.simplify_go <- function(x, similarity_cutoff = 0.9){

  go_enrichment_result_only_enriched <-
    x %>%
    purrr::map(~ clusterProfiler::simplify(..1, cutoff = similarity_cutoff))

  return(go_enrichment_result_only_enriched)
}


#' Generate a list of GO EMAP plot.
#'
#' @param x a list containing output of GO enrichment. Each element of the list is an output of  [clusterProfiler::enrichGO()].
#' @param n_terms number of terms to show.
#' @param color_terms_by
#'
#' @return
#' @export
#'
#' @keywords internal
.generate_go_emap_plot <- function(x, n_terms = 30 , color_terms_by = "p.adjust"){

  safely_plot <- purrr::safely(.f = function(x, y,...){
    ep <- enrichplot::emapplot(x,...)
    ep <-   ep + ggplot2::ggtitle(y)
    return(ep)
  } , otherwise = NULL)

  go_plots <- purrr::map(names(x) , ~
                           safely_plot(x = x[[..1]],
                                       y = ..1,
                                       showCategory = n_terms,
                                       color = color_terms_by,
                                       layout = "kk"))

  # simplify list
  go_plots <- go_plots %>%
    purrr::transpose() %>%
    purrr::simplify_all()

  # keep result
  go_plots <- go_plots$result

  # assign names
  names(go_plots) <- names(x)

  # remove element if NULL
  go_plots <- go_plots %>%
    purrr::discard(is.null)

  go_plots

}


#' Group replicates by sample - returns list
#' @rdname get_replicates_by_sample
#'
#' @keywords internal
#' @export
.get_replicates_by_sample_list <- function(x){

  # validate x
  .validate_parcutils_obj(x)

  y <- .group_replicates_by_sample(x)

  y_grpd <-
    y %>%
    dplyr::group_by(groups) %>%
    tidyr::nest() %>%
    dplyr::mutate(data = purrr::map(data, ~ ..1 %>% dplyr::pull(1)))

  y <-
    y_grpd %>%
    dplyr::pull(2)

  names(y) <- y_grpd %>% dplyr::pull(1)

  return(y)
}



#' Map metadata (GC, length and seq) to GRanges object
#'
#' @param x an object of class GRanges
#' @param bs_genome_object an object of class BSgenome, default \code{BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38}
#'
#' @return
#' @export
#'
#' @keywords internal
.map_granges_metadata <- function(x, bs_genome_object = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38){

  stopifnot("x must be the object of class GRanges" = is(x, "GRanges"))
  stopifnot("bs_genome_object must be an object of class BSgenome" = is(bs_genome_object, "BSgenome"))

  x <- x %>%
    # add sequence for each intron
    dplyr::mutate(seq =  BSgenome::getSeq(bs_genome_object,x)) %>%

    # add GC for each intron
    dplyr::mutate(GC =  BSgenome::letterFrequency(seq, letters = "GC", as.prob = T) %>%  as.numeric()) %>%

    # add length for each intron
    dplyr::mutate(length =  BSgenome::width(seq))

  return(x)
}


#' For each intron in the object irFinderSdata assign unique intron id
#'
#' @param x an object of class irFinderSdata
#'
#' @return an object of class irFinderSdata
#' @export
#' @keywords internal
#' @examples
#'
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' .parcutils_assign_intron_identifier(x)
.parcutils_assign_intron_identifier  <- function(x){
  # assign intron id to each element in the x

  .validate_irfinders_object(x)
  x <- purrr::map(x  , ~ ..1 %>% dplyr::mutate(intron_id = stringr::str_c("intron", 1:dplyr::n(), sep = "_")))
  x <- .assign_class_irFinderSdata(x)
  return(x)
}

#' Check if the object belongs to class irFinderSdata.
#'
#' @param x an object to be tested for class 'irFinderSdata'
#'
#' @return error.
#' @export
#'
#' @keywords internal
.validate_irfinders_object <- function(x){
  stopifnot("x must be an object of class irFinderSdata" = is(x, "irFinderSdata") )
}


#' Assign class irFinderSdata
#' @description This function does all mandatory checks before it assigns class irFinderSdata
#' @param x a list or dataframe to which class irFinderSdata to assign.
#'
#' @return an object of class irFinderSdata
#' @export
#' @keywords internal
.assign_class_irFinderSdata <- function(x){

  # x can be a dataframe or list of dataframes
  # if dataframe it must have mandatory columns
  # if a list it must have mandatory columns in each dataframe and same number of rows in each dataframe

  if(is(x , "data.frame")){
    .check_mendate_columns_for_irFinderSdata(x)
    x <- list(x)
    class(x) <- c("irFinderSdata", class(x))
  }

  if(is(x , "list")){
    purrr::walk(x, ~ .check_mendate_columns_for_irFinderSdata(..1))

    # all elems of the list have same number of rows
    x_nrows <- purrr::map_int(x, ~..1 %>% nrow())

    if(!all(x_nrows==x_nrows[1])){
      stop("All elements in x must have same number of rows.")
    }

    class(x) <- c("irFinderSdata", class(x))
  }

  return(x)

}



#' Check mandatory columns in a dataframe
#' @description This function helps to create mandatory column in a dataframe for the object irFinderSdata.
#' @param x dataframe
#' @param mandat_columns a character vector denoting mandatory columns in x.
#'
#' @return TRUE or ERROR
#' @export
#' @keywords internal
.check_mendate_columns_for_irFinderSdata <- function(x,
                                                     mandat_columns = c("chr",
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

  stopifnot("x must be a dataframe" = is(x , "data.frame"))

  x_cols <- colnames(x)

  col_not_found <- mandat_columns[is.na(match(mandat_columns, x_cols))]
  if(length(col_not_found) > 0){
    stop(cli::format_error(c("x" ="Column{?s} {col_not_found} not found")))
  }

  return(TRUE)

}


#' Get intron annotations
#'
#' @param x an object of class irFinderSdata
#' @param intron_ids a character vector of intron ids for which
#' annotations to be retrieved.
#'
#' @return a dataframe with columns \code{chr, start, end, intron_id,
#' null, strand, gene_id, gene_symbol}
#' @export
#' @keywords internal
#' @examples
#' example_dir <- system.file("extdata" ,"ir_data" , package="parcutils")
#' example_files <- fs::dir_ls(example_dir, glob = "*-IR-nondir*.txt")
#' names(example_files) <- stringr::str_replace(example_files , pattern = ".*/", replacement = "") %>%
#' stringr::str_replace(pattern = "-IR-nondir.txt", replacement = "")
#' x <- read_irfinderS_output(files = example_files,  add_prefix_chr = FALSE)
#' .get_intron_annotations(x, intron_ids = paste("intron", sample(1:100, 10),sep ="_"))
.get_intron_annotations <- function(x, intron_ids){

  .validate_irfinders_object(x)

  stopifnot("'intron_ids' must be a character vector." = is.character(intron_ids))

  # first six columns of irfinders output
  dd <- x[[1]]  %>%
    dplyr::select(1:6 , intron_id)  %>%
    dplyr::filter( intron_id %in% intron_ids) %>%
    dplyr::mutate(gene_id = stringr::str_replace(name, pattern = "(.*)/(.*)/(.*)",replacement = "\\2")) %>%
    dplyr::mutate(gene_symbol = stringr::str_replace(name, pattern = "(.*)/(.*)/(.*)",replacement = "\\1")) %>%
    dplyr::mutate(name = intron_id) %>%
    dplyr::select(-intron_id) %>%
    dplyr::rename(intron_id = name) %>%
    dplyr::distinct()

  return(dd)
}



#' Prepare object of parcutils_ir from the object of parcutils
#'
#' @param x an object of class 'parcutils'.
#'
#' @return
#' @export
#' @keywords internal
.prepare_parcutils_ir <-  function(x){
  .validate_parcutils_obj(x)
  y <- x %>% tibble::tibble()
  id_column <- y$norm_counts[[1]][[1]] %>% colnames() %>% .[1]
  # fix 1st column names

  y <- y %>%
    dplyr::mutate(norm_counts = purrr::map(norm_counts, function(.x){
      purrr::map(.x, ~ ..1 %>% dplyr::rename( !!id_column := 1))})) %>%
    dplyr::mutate(dsr_tibble = purrr::map(dsr_tibble, function(.x){
      .x %>% dplyr::rename(!!id_column := 1)}))  %>%
    dplyr::mutate(dsr_tibble_deg = purrr::map(dsr_tibble_deg, function(.x){
      .x %>% dplyr::rename(!!id_column := 1)
    }))

  class(y) <- c("parcutils_ir",class(x))
  return(y)
}



#' Save DEG results in an excel file.
#'
#' @param x an object of class parcutils.
#' @param file a file path.
#'
#' @return
#' @export
#' @keywords internal
.save_deg_results <- function(x, file){

  .validate_parcutils_obj(x)
  de_table <- x$dsr_tibble_deg
  ret <- writexl::write_xlsx(de_table,path =  file)
  return(ret)
}

#' Save DE ir results in an excel file.
#'
#' @param x an object of class parcutils_ir
#' @param file a file path.
#'
#' @return a file path.
#' @export
#' @keywords internal
.save_de_ir_results <- function(x , file ){

  .validate_parcutils_ir_object(x)

  oo <- purrr::map(x$comp , ~{
    de_table <- x$dsr_tibble_deg[[.x]]
    intron_annot_table <-  x$intron_annot[[.x]]

    # join de and annotation tables.
    intron_annot_table %>%
      dplyr::left_join( de_table ,by = "intron_id") %>%
      # add IGV search text
      dplyr::mutate(igv_search_text = stringr::str_c(chr, ":", start, "-",end))
  })
  names(oo) <- x$comp
  ret <- writexl::write_xlsx(oo,path =  file)
  return(ret)
}

#' Validate parcutils_ir object
#'
#' @param x
#'
#' @return
#' @export
#' @keywords internal
.validate_parcutils_ir_object <- function(x) {
  stopifnot("x must be an object of class 'parcutils_ir'.
           Usually x is derived by parcutils::run_deseq_analysis_ir()." = is(x, "parcutils_ir"))
}

