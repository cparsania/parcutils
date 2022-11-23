
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
#' \dontrun{
#' }
.piarwise_upset <- function(x, sample_comparison , color_up = "#b30000", color_down = "#006d2c",... ){

  # validate sample_comparison

  stopifnot("sample_comparison must be a character vector of length 2." = is.character(sample_comparison) & length(sample_comparison) == 2)

  yy <- purrr::map(sample_comparison %>%
                     purrr::set_names(sample_comparison), ~parcutils::get_genes_by_regulation(x = x, sample_comparison = ..1,regulation = "both") ) %>%
    purrr::map(~split(names(..1), ..1))

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
                           keep.order = TRUE,
                           text.scale = 1.5,
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
#' @param var_gene_id // TO DO
#' @param var_feature_type // TO DO
#' @param var_feature_biotype // TO DO
#' @param var_gene_name // TO DO, to be implement yet.
#'
#' @return // TO DO
#' @keywords internal
.filter_gff <- function(gtf_file ,

                       feature_type = "gene" ,
                       feature_biotype = "protein_coding",

                       var_gene_name = "gene_name",
                       var_gene_id = "gene_id",
                       var_feature_type = "type" ,
                       var_feature_biotype = "gene_biotype"){

  # gtf_file = gff
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
                                                  sample_comparisons = x$de_comparisons[1])

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
#' @return a list of GO plots.
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
#' @return a dataframe.
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
    dplyr::mutate(GC =  BSgenome::letterFrequency(seq, letters = "GC", as.prob = TRUE) %>%  as.numeric()) %>%

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
#' @return an object of class 'parcutils_ir'
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
#' @return output of \code{writexl::write_xlsx()}
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

  oo <- purrr::map(x$de_comparisons , ~{
    de_table <- x$dsr_tibble_deg[[.x]]
    intron_annot_table <-  x$intron_annot[[.x]]

    # join de and annotation tables.
    intron_annot_table %>%
      dplyr::left_join( de_table ,by = "intron_id") %>%
      # add IGV search text
      dplyr::mutate(igv_search_text = stringr::str_c(chr, ":", start, "-",end))
  })
  names(oo) <- x$de_comparisons
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


#' Prepare count data for the function [parcutils::run_deseq_analysis()]
#'
#' @param counts a character sting or dataframe.
#'
#' @return a dataframe.
#' @export
#' @keywords internal
.get_count_data  <- function(counts){
  if(!inherits(counts , what = c("character","tbl_df","tbl" ,"data.frame"))){
    cli::cli_abort("{.arg counts} must be an object of class {.cls character} or {.cls data.frame}")
  }

  if(inherits(counts , what = c("character"))){
    cli::cli_alert_info(" Reading count file ...")
    count_data <- readr::read_delim(counts,
                                    delim = delim,
                                    comment = comment_char,
                                    trim_ws = TRUE,
                                    col_names = TRUE,
                                    show_col_types = FALSE ,
                                    progress = TRUE)
    cli::cli_alert_success("Done.")
  } else {
    count_data <- counts
  }
  # remove duplicate rows if any
  count_data %<>% dplyr::distinct()

  return(count_data)
}


#' Prepare sample information for the function [parcutils::run_deseq_analysis()]
#'
#' @param sample_info a character string or data frame.
#'
#' @return a data frame
#' @export
#'
#' @keywords internal
.get_sample_information <- function(sample_info){

  # value for sample_info either from a class character, tbl_df, tbl or data.frame

  if(!inherits(sample_info, what = c("character","tbl_df","tbl" ,"data.frame"))){
    cli::cli_abort("{.arg sample_info} must be an object of class {.cls character} or {.cls data.frame}")

  }

  # if sample_info is a character string then read data from file.

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
      cli::cli_abort("{.arg sample_info} must contains two columns deliminated by '\\t'")
    }
  } else{
    sample_info_data = sample_info_data
  }

  # remove duplicates
  sample_info_data %<>% dplyr::distinct()

  return(sample_info_data)

}



#' Convert NA to 0 from into count data for the function [parcutils::run_deseq_analysis()]
#'
#' @param count_data a dataframe
#'
#' @return a dataframe
#' @export
#'
#' @keywords internal
.count_data_convert_NA_to_0 <- function(count_data){
  count_data_select_NA <- count_data %>%
    TidyWrappers::tbl_keep_rows_NA_any()

  if(nrow(count_data_select_NA) >= 1) {
    cli::cli_alert_warning("Following gene(s) contains NA in one or more samples. Either remove them or Default NA will be converted to 0.")
    genes_NA <- count_data_select_NA %>%
      dplyr::pull(!!column_geneid_quo)
    cli::cli_alert(genes_NA)
    # replace NA with 0
    count_data %<>% dplyr::mutate_if(is.numeric , ~ ..1 %>% tidyr::replace_na(0))
  }
  return(count_data)
}



#' Filter count data using user thresholds for the function [parcutils::run_deseq_analysis()]
#'
#' @param count_data a data frame
#' @param column_geneid a character string
#' @param sample_info a data frame
#' @param min_counts a numeric value
#' @param min_replicates a numeric value
#'
#' @return
#' @export
#'
#' @keywords internal
.filter_rows_from_count_data <- function(count_data, column_geneid, sample_info,
                                         min_counts,min_replicates){

  # For DESeq analysis genes having non zero count in at least one sample will be used.
  # In addition genes will be filtered based on user defined cutoffs -
  # 1. minimum counts to each gene
  # 2. minimum number of replicates with minimum counts

  # quote column names
  column_geneid_quo <- rlang::enquo(column_geneid)

  sample_info_colnames <- colnames(sample_info)

  genes_to_keep <- count_data %>%

    # convert data wide to long format.
    tidyr::gather(-!!column_geneid_quo, key = "sample", value = "counts") %>%

    # join sampleinfo table
    dplyr::left_join(sample_info, by = c("sample" = sample_info_colnames[1])) %>%

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

  count_data_filter <- count_data %>%
    dplyr::filter(!!rlang::sym(column_geneid) %in% genes_to_keep)

  return(count_data_filter)
}


#' Check availability of sample comparisons in the object of parcutils
#'
#' @param x an object of class \code{parcutils}
#' @param sample_comparisons a character vector
#'
#' @return either TRUE or abort.
#' @export
#'
#' @keywords internal
.is_sample_comparsions_present_in_obj_parcutils <- function(x, sample_comparisons){

  .validate_parcutils_obj(x)
  not_present <- sample_comparisons[!(sample_comparisons %in% x$de_comparisons)]
  if(length(not_present) > 0 ){
    cli::cli_abort("{.arg sample_comparisons} - {.emph {cli::col_red({not_present})}} {?is/are} not present in {.arg x}")
  } else{
    return(TRUE)
  }
}


#' Check availability of genes in the object of parcutils
#'
#' @param x an object of class \code{parcutils}.
#' @param genes a character vector.
#'
#' @return either TRUE or abort.
#' @export
#'
#' @keywords internal
.is_genes_present_in_obj_parcutils <- function(x, genes){
  all_genes <- .get_all_expressed_genes(x) %>% names()

  # check if values from labels present in x
  not_present <- genes[!genes %in% all_genes]
  if(length(not_present) >0 ){
    cli::cli_abort("{.arg genes} - {.emph {cli::col_red({not_present})}} {?is/are} not present in {.arg x}")
  } else{
    return(TRUE)
  }

}


#' Get gene groups for fc scatter plot
#'
#' @param x an object of class parcutils
#' @param sample_comparisons a character vector
#'
#' @return a dataframe
#' @export
#'
#' @keywords internal
.get_gene_groups_for_fc_scatter_plot <- function(x, sample_comparisons,color_label){

  # prepare groups to color data points
  genes_by_regul <- parcutils::get_genes_by_regulation(x = x,
                                                       sample_comparisons =sample_comparisons,
                                                       regulation = "all",simplify = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$regul,.data$genes) %>%
    dplyr::add_count(name = "counts") %>%
    dplyr::ungroup()

  # assign groups -> both_up, both_down, both

  if(color_label == "both_up"){
    gene_groups <- genes_by_regul %>%
      dplyr::mutate(group = dplyr::case_when(counts == 2 & grepl("up", regul) ~
                                              stringr::str_c("both_",regul,sep = ""),
                                            TRUE ~ "other")) %>%
      dplyr::select(.data$genes, .data$group) %>% dplyr::distinct()

  } else if(color_label == "both_down"){
    gene_groups <- genes_by_regul %>%
      dplyr::mutate(group = dplyr::case_when(counts == 2 & grepl("down", regul) ~
                                              stringr::str_c("both_",regul,sep = ""),
                                            TRUE ~ "other")) %>%
      dplyr::select(.data$genes, .data$group) %>% dplyr::distinct()
  } else {
    gene_groups <- genes_by_regul %>%
      dplyr::mutate(group = dplyr::case_when(counts == 2 & grepl("up|down", regul) ~
                                              stringr::str_c("both_",regul,sep = ""),
                                            TRUE ~ "other")) %>%
      dplyr::select(.data$genes, .data$group) %>% dplyr::distinct()
  }

}


#' Perform k-means clustering
#'
#' @param data a dataframe. First column will be used as row identifier.
#' Therefore, values in the first column must be unique.
#' @param km an integer denoting number of clusters to return.
#'
#' @return a dataframe with two columns -1) 1st column from `data` and 2) clusters
#' @export
#'
#' @keywords internal
.perform_kemans_cluster <- function(data, km ){

  data_column_names <- data %>% colnames()
  data_id_col <- data_column_names[1]

  # convert data to matrix
  data_mat <- data %>%
    as.data.frame() %>%
    tibble::column_to_rownames(data_id_col) %>% as.matrix()
  set.seed(123)
  km_out <- kmeans(data_mat, centers = km )
  ret <- km_out$cluster %>%
    tibble::enframe() %>%
    dplyr::rename(!!rlang::sym(data_id_col) := name, clusters = value) %>%
    dplyr::mutate(clusters = stringr::str_c("clust_",clusters, sep = ""))
  return(ret)
}


#' Generate a gene expression line plot.
#' @description This function called internally from the function
#' `get_gene_expression_line_plot`.
#' @param line_plot_data_wide a dataframe containing data for line plot.
#' @param km a numeric, indicates number of clusters for k-means clustering.
#' @param facet_clusters a logical, denotes whether to facet clusters generated by k-means.
#' @param scale_log10 a logical, denotes whether to transform Y axis on log10 scale.
#' @param line_transparency a numeric, denotes value for line transparency.
#' @param show_average_line a logical, denotes whether to show average line.
#' @param average_line_color a character string, denotes color for average line.
#' @param average_line_size a numeric, denotes size for average line.
#' @param average_line_summary_method character string either \code{mean} or  \code{median},
#' denotes summary method for average line.
#'
#' @return a ggplot2.
#' @export
#'
#' @keywords internal
.generate_a_lineplot <- function(line_plot_data_wide,
                                                  km,
                                                  facet_clusters,
                                                  scale_log10,
                                                  line_transparency,
                                                  show_average_line,
                                                  average_line_color,
                                                  average_line_size,
                                                  average_line_summary_method
){

  gene_id_col <- line_plot_data_wide %>% colnames() %>%.[1]

  # long data
  line_plot_data_long <- line_plot_data_wide %>%
    tidyr::pivot_longer(-!!rlang::sym(gene_id_col), names_to = "samples", values_to = "values") %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(row_num =  1:dplyr::n())

  # perform clustering
  if(!is.null(km)){
    km_output <- .perform_kemans_cluster(data = line_plot_data_wide, km = km)
    line_plot_data_long <- line_plot_data_long %>%
      dplyr::left_join(km_output, by = gene_id_col)
  }

  # plot
  gp <- line_plot_data_long %>%
    ggplot2::ggplot(ggplot2::aes(x = samples, y = values)) +
    ggplot2::geom_point()

  # add layers
  if(is.null(km)){
    gp <- gp + ggplot2::geom_line(alpha = line_transparency, ggplot2::aes(group = row_num))
  } else {
    gp <- gp + ggplot2::geom_line(alpha = line_transparency, ggplot2::aes(group = row_num, col = clusters))
  }

  # modify theme and scale
  if(TRUE){
    gp <- gp +   ggplot2::theme(text = ggplot2::element_text(size = 12)) +
      ggplot2::theme_bw() +
      ggplot2::scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_long_scale()))
  }
  # scale log10
  if(scale_log10) {
    suppressMessages(
      gp <- gp +
        ggplot2::scale_y_log10(labels = scales::label_number(scale_cut = scales::cut_long_scale()))
    )
  }
  # change label
  if(TRUE){
    gp  <- gp + ggplot2::ylab("Normalised Gene Expression") +
      ggplot2::xlab("Samples")
  }

  # show average line
  if(show_average_line){
    if(!is.null(km)){
      avg_line_data <- line_plot_data_long %>%
        dplyr::group_by(samples, clusters) %>%
        dplyr::summarise_at(.vars = dplyr::vars(values),
                            .funs = average_line_summary_method)
      gp <- gp +
        ggplot2::geom_line(data = avg_line_data %>%
                             dplyr::mutate(row_num = dplyr::row_number()) ,
                           ggplot2::aes( group = row_num),
                           col = average_line_color ,
                           size = average_line_size)

    } else {
      avg_line_data <- line_plot_data_long  %>%
        dplyr::group_by(samples) %>%
        dplyr::summarise_at(.vars = dplyr::vars(values),
                            .funs = average_line_summary_method)

      gp <- gp +
        ggplot2::geom_line(data = avg_line_data %>%
                             dplyr::mutate(row_num = dplyr::row_number()) ,
                           ggplot2::aes( group = 1),
                           col = average_line_color ,
                           size = average_line_size)
    }

  }
  # facet
  if(facet_clusters){
    if(!is.null(km)){
      gp <- gp + ggplot2::facet_wrap(~clusters)
    }
  }
  return(gp)
}




#' Get the gene names stored in the parcutils object
#'
#' @param x an object of class parcutils
#' @param query_genes a character vector denoting genes for which corresponding gene names to be derived from the parcutils object.
#'
#' @return a dataframe containing two columns - 1) queried genes and 2) gene names from the object parcutils.
#' @export
#'
#' @examples
#'
#' res = .get_parcutils_object_example()
#' .get_parcutils_obj_gene_names(x = res, query_genes= c("LSS", "GOLGA8S", "FXYD5", "CDC34"))
#'
#' @keywords internal
.get_parcutils_obj_gene_names <- function(x, query_genes){


  lookup_table <- parcutils:::.get_all_expressed_genes(x = x) %>%
    names()

  names(query_genes) <- query_genes

  purrr::map(query_genes , ~which(stringr::str_detect(pattern = ..1,string = lookup_table))) %>%
    tibble::enframe(value = "index") %>%
    tidyr::unnest(keep_empty = T) %>%
    dplyr::mutate(full_name = lookup_table[index]) %>%
    dplyr::select(-index)


}


#' Obtain an example object of the class parcutils
#'
#' @return an object of class \code{parcutils}
#' @export
#'
#' @keywords internal
.get_parcutils_object_example <- function(){
  count_file <- system.file("extdata","toy_counts.txt" , package = "parcutils")
  count_data <- readr::read_delim(count_file, delim = "\t",show_col_types = FALSE)

  sample_info <- count_data %>% colnames() %>% .[-1]  %>%
    tibble::tibble(samples = . , groups = rep(c("control" ,"treatment1" , "treatment2"), each = 3) )


  res <- run_deseq_analysis(counts = count_data ,
                            sample_info = sample_info,
                            column_geneid = "gene_id" ,
                            cutoff_lfc = 1 ,
                            cutoff_pval = 0.05,
                            group_numerator = c("treatment1", "treatment2") ,
                            group_denominator = c("control"))

  return(res)
}


