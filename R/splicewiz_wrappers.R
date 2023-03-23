#' Generate a barplot showing counts of alternate splice events from the object NxtSe.
#' @description This plot is useful to show counts of filtered spliced events once spliceWiz filters are applied.
#' @param x an object of class \code{NxtSE}. Usually \code{se} or \code{se.filtered}.
#' @param elem_text_size an integer denoting size of the various text elements in the plot.
#' @param text_count_size an integer denoting size of the counts mentioned above each bar.
#'
#' @return ggplot.
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' get_ASE_counts_barplot(se, elem_text_size = 15, text_count_size = 4)
#'
#'
#' @keywords internal
get_ASE_counts_barplot <- function(x, elem_text_size = 15, text_count_size = 10){

  stopifnot("x must be the class of NxtSE" = is(x, "NxtSE"))
  gp <- x@elementMetadata %>%
    tibble::as_tibble() %>%
    dplyr::select(EventName, EventType) %>%
    dplyr::group_by(EventType) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(EventType = forcats::fct_reorder(EventType, -n)) %>%
    ggplot(aes(x = EventType, y = n)) +
    geom_bar(position ="dodge", stat = 'identity') +
    geom_text(aes(label=n), inherit.aes = TRUE, position = position_dodge(width = 0.9),size = text_count_size, vjust=-0.25) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = elem_text_size)) +
    ggplot2::ylab("Counts") +
    ggplot2::xlab("EventType") +
    ggeasy::easy_rotate_labels(which = "x")

  return(gp)

}



#' Perform ASE diff analysis using edgeR.
#'
#' @param x an object of class \code{NxtSe}.
#' @param test_factor refer to the argument \code{test_factor} in [SpliceWiz::ASE_edgeR()].
#' @param test_nom,test_denom A pair of vectors, usually the same length, denoting nominator and denominator conditions to test for differential ASE, respectively.
#' Usually the "treatment" condition is used in nominator and "control" is used in the denominator.
#' @param regul_based_upon one of the numeric choices  1, 2, or 3. Default 1 i.e. categorized diff. ASE by pvalue and log2fc.
#' ## if 1 ...
#'  - Up : log2fc >= cutoff_lfc & pvalue <= cutoff_pval
#'  - Down : log2fc  <= (-1) * cutoff_lfc & pvalue <= cutoff_pval
#'  - Other : remaining genes
#' ## if 2 ...
#'  - Up : log2fc >= cutoff_lfc & padj <= cutoff_padj
#'  - Down : log2fc  <= (-1) * cutoff_lfc & padj <= cutoff_padj
#'  - Other : remaining genes
#' ## if 3 ...
#'  - Up : log2fc >= cutoff_lfc & pvalue <= cutoff_pval & padj <= cutoff_padj
#'  - Down : log2fc  <= (-1) * cutoff_lfc & pvalue <= cutoff_pval & padj <= cutoff_padj
#'  - Other : remaining genes
#' @param cutoff_lfc minimal threshold for log2fold change, default 1 (2 fold).
#' @param cutoff_pval minimal threshold for pvalue, default 0.05. P-value threshold will be applied only when
#'  `regul_based_upon` is either 1 or 3.
#' @param cutoff_padj minimal threshold for Padj, default 0.01. Padj threshold will be applied only when
#'  `regul_based_upon` is either 2 or 3.
#' @param ... Other parameters passed to [SpliceWiz::ASE_edgeR()].
#'
#' @return an object of class parcutils_se.
#' @export
#' @importFrom rlang inject
#' @importFrom purrr map2
#' @importFrom SpliceWiz ASE_edgeR
#' @importFrom stringr str_c
#' @importFrom rlang `!!!`
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated")
#' run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_padj = 1, regul_based_upon = 2)
#'
run_ase_diff_analysis <- function(x, test_nom ,test_denom, test_factor = "condition",cutoff_lfc = 1, cutoff_pval = 0.05 , cutoff_padj= 0.01, regul_based_upon = 1, ...){

  # get all ... arguments
  args = c(...)

  # run edger.
  res_ase_diff_raw  <- purrr::map2(test_nom, test_denom , ~ {
    rlang::inject(SpliceWiz::ASE_edgeR(se=x,
                                       test_factor = test_factor,
                                       test_nom = ..1 ,
                                       test_denom = ..2, !!! args))
  }
  )

  # assign names to each comparison
  names(res_ase_diff_raw) <- stringr::str_c(test_nom, test_denom,sep = "_VS_")

  # convert edger result into the tibble and fix column names

  res_ase_diff_tibble <- purrr::map(names(res_ase_diff_raw), ~ res_ase_diff_raw[[..1]] %>% tibble::as_tibble())
  names(res_ase_diff_tibble) <- names(res_ase_diff_raw)

  res_ase_diff_tibble <- purrr::map(res_ase_diff_tibble , ~.fix_splicewiz_res_edger_columns(..1) %>%
    dplyr::select(c("event_name","event_type","event_region"), dplyr::contains("avg_psi"), log2FoldChange,pvalue, padj))

  # get matrix of AVG_PSI for each sample into the comparison

  avg_psi <- purrr::map(res_ase_diff_raw , ~..1 %>% .fix_splicewiz_res_edger_columns()  %>% dplyr::select(event_name, event_type, event_region, dplyr::contains("avg_psi")))

  # mark diff regulation "Up" "Down" and "Other" in the res_ase_diff_tibble

  res_ase_diff_tibble_annot <- purrr::map(res_ase_diff_tibble, ~ .categorize_diff_ase(..1,log2fc_cutoff = cutoff_lfc, pval_cutoff = cutoff_pval, padj_cutoff = cutoff_padj, regul_based_upon = regul_based_upon))

  # get diff ase counts

  res_ase_diff_summary <-   purrr::map(res_ase_diff_tibble_annot , ~ ..1 %>%
                                         dplyr::group_by(regul) %>%
                                         dplyr::tally() ,.id = "cond")

  # Combine everything into a tibble and identify it as an object of class "parcutils_ase".

  ret_obj <- tibble::tibble(de_comparisons = names(res_ase_diff_raw)) %>%
    dplyr::mutate(numerator = stringr::str_replace(de_comparisons, pattern = "_VS.*", replacement = "")) %>%
    dplyr::mutate(denominator = stringr::str_replace(de_comparisons, pattern = ".*VS_", replacement = "")) %>%
    dplyr::mutate(avg_psi = avg_psi) %>%
    dplyr::mutate(res_ase_diff_raw = res_ase_diff_raw) %>%
    dplyr::mutate(res_ase_diff_tibble = res_ase_diff_tibble) %>%
    dplyr::mutate(res_ase_diff_tibble_annot = res_ase_diff_tibble_annot) %>%
    dplyr::mutate(res_ase_diff_summary = res_ase_diff_summary)


  class(ret_obj) <- c("parcutils_ase" , class(ret_obj))

  return(ret_obj)

}


#' Generate a barplot of diff ASE counts.
#'
#' @param x an object of class parcutils_ase.
#' @param col_up a character string, default `#a40000`, denoting valid a color name for "Up" regulated genes.
#' @param col_down a character string, default `#16317d`, denoting valid a color name "Down" regulated genes.
#' @param font_size a numeric, default 12, denoting font size in the plot.
#' @param show_counts a logical, default `TRUE`, denoting whether to show counts on each bar.
#'
#' @return a bar plot.
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_padj = 1, regul_based_upon = 2)
#' get_diff_ASE_count_barplot(res)
get_diff_ASE_count_barplot <- function(x,
                                       col_up="#a40000",
                                       col_down="#16317d",
                                       font_size = 12,
                                       show_counts = TRUE){

  .validate_parcutils_ase_obj(x)

  gp <- x$res_ase_diff_summary %>%
    tibble::enframe(name = "comparison" ,
                    value = "deg_count") %>%
    tidyr::unnest(cols = "deg_count") %>%
    dplyr::filter(regul != "other") %>%
    dplyr::mutate(regul = forcats::fct_relevel(regul, c("Up","Down"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = regul, y = n , fill = regul)) +
    ggplot2::facet_wrap(~comparison) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(breaks = c("Up","Down"),
                               values = c(col_up,col_down) ) +
    ggplot2::theme(text = ggplot2::element_text(size = font_size)) +
    ggplot2::xlab("Regulation") +
    ggplot2::ylab("Counts") +
    ggplot2::guides(fill =  ggplot2::guide_legend("Regulation"))

  if(show_counts){
    gp <- gp + ggplot2::geom_text(ggplot2::aes(label = n),
                                  inherit.aes = T,
                                  position = position_dodge(width = 0.9),
                                  size = font_size/3, vjust = 0.05, col = "black")
  }

  return(gp)
}


#' Get a data matrix (PSI, logit, z-score) for a given ASE
#'
#' @param se an object of class NxtSE.
#' @param event_names a character vector denoting valid ASE names to plot in the heatmap.
#' @param samples a character vector denoting valid sample names to plot in the heatmap.
#' @param replicates logical, default TRUE, denoting whether to show sample replicates or not.
#' @param method a character string denoting a method for data value. Can be one of the "PSI", "logit" or "Z-Score"
#' @param column_condtion a name of column in se storing condition. Default: "condition"
#'
#' @return a data matrix
#' @export
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_padj = 1, regul_based_upon = 2)
#' event_names = get_ASEsets_by_regulation(x = res, sample_comparisons = "A_VS_B", regul = "all") %>% unlist()
#' get_ASE_data_matrix(se, event_names , samples = c("A", "B"), column_condition ="treatment")
#'
get_ASE_data_matrix <- function(se , event_names,samples,replicates=TRUE, method ="PSI", column_condition ="condition"){

  samples_repli <- SpliceWiz::colData(se) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "samples") %>%
    dplyr::filter(!!rlang::sym(column_condition) %in% !!samples) %>%
    dplyr::pull(samples)

  data_matrix <- SpliceWiz::makeMatrix(se = se,
                                       event_list = event_names,
                                       sample_list = samples_repli,
                                       method = method) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "event_name") %>%
    tibble::tibble()


  return(data_matrix)
}



#' Generate a heatmap for a given list of ASE.
#'
#' @param se an object of class NxtSE.
#' @param event_names a character vector denoting valid spliceWiz event names.
#' @param samples a character vector denoting valid sample names.
#' @param column_condition a name of column in se storing condition. Deafault: "condition".
#' @param method one of the character strings: "Z-score", "PSI" ,"logit".
#' @param show_replicates logical, default TRUE, denoting whether to show sample replicates or not.
#' @param show_row_names logical, whether to show row names in the heatmap.
#' @param cluster_rows logical, whether to clusters rows in the heatmap.
#' @param cluster_columns logical, whether to cluster column in the heatmap or not.
#' @param show_column_dend logical, whether to show column dendrogram in the heatmap or not.
#' @param show_row_dend logical, whether to show row dendrogram in the heatmap or not.
#' @param ... other parameters to pass [ComplexHeatmap::Heatmap()]
#'
#' @return a heatmap
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_padj = 1, regul_based_upon = 2)
#' event_names = get_ASEsets_by_regulation(x = res, sample_comparisons = "A_VS_B", regul = "all") %>% unlist()
#'
#' get_ase_data_matrix_heatmap(se, event_names = event_names, samples = c("A" ,"B"), column_condition    = "treatment")
#' get_ase_data_matrix_heatmap(se, event_names = event_names, samples = c("A" ,"B"), column_condition    = "treatment",method = "Z-score", cluster_rows = TRUE)
get_ase_data_matrix_heatmap <- function(se,
                                 event_names,
                                 samples,
                                 column_condition = "condition",
                                 method = "PSI",
                                 show_replicates=TRUE, # not implemented yet
                                 show_row_names = FALSE,
                                 cluster_rows = FALSE,
                                 cluster_columns = FALSE,
                                 show_column_dend = FALSE,
                                 show_row_dend = FALSE,...){


  dat <- get_ASE_data_matrix(se = se,
                               event_names = event_names ,
                               column_condition = column_condition,
                               samples = samples,method = method)

  dat <- dat %>%
    TidyWrappers::tbl_remove_rows_NA_any() %>% as.data.frame() %>%
    tibble::column_to_rownames("event_name")

  ComplexHeatmap::Heatmap(dat, show_row_names = show_row_names,
                          cluster_rows = cluster_rows,
                          show_column_dend = show_column_dend,
                          cluster_columns = cluster_columns,
                          show_row_dend = show_row_dend, ...)


}



#' Get gene names from spliceWiz's EventName column.
#'
#' @param x a character vector denoting EventNames of SpliceWiz.
#'
#' @return a character vector of gene names
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' res_edgeR <- SpliceWiz::ASE_edgeR(se, "treatment", "A", "B", useQL = FALSE)
#' get_genes_from_event_name(res_edgeR$EventName)
get_genes_from_event_name <- function(x){

  x %>% stringr::str_replace( pattern = ".*;(.*)-.*-.*", replacement = "\\1") %>%
    stringr::str_replace(pattern = "\\/.*", replacement = "")
}


#' Convert SpliceWiz's EventRegion column to GRanges object.
#'
#' @param event_region a character vector denoting eventregion from spliceWiz out. Typically looks like this: "14:65968689-65971731/+"
#' @param prefix a character string denoting a prefix to be added to chromosome names. Default none. Typically the string can be "chr" or "Chr"
#'
#' @return an object of class [GenomicRanges::GRanges()]
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' event_region <- se@elementMetadata$Event1a %>% sample(10)
#' event_region_to_granges(event_region = event_region)
#'
event_region_to_granges <- function(event_region, prefix = ""){

  col_seqnames = glue::glue("seqnames")
  col_start = glue::glue("start")
  col_end = glue::glue("end")
  col_strand = glue::glue("strand")

  tbl <- tibble::tibble(event_region = event_region) %>%
    dplyr::mutate(!!col_seqnames := stringr::str_replace(event_region, ":.*","")) %>%
    dplyr::mutate(!!col_seqnames := stringr::str_c(prefix,!!rlang::sym(col_seqnames),sep = "")) %>%
    dplyr::mutate(!!col_start:= stringr::str_replace(event_region, ".*:(\\d+)-.*","\\1") %>%
                    as.numeric()) %>%
    dplyr::mutate(!!col_end := stringr::str_replace(event_region, ".*:(\\d+)-(\\d+)\\/.*","\\2") %>% as.numeric()) %>%
    dplyr::mutate(!!col_strand := stringr::str_replace(event_region, ".*:(\\d+)-(\\d+)\\/(.*)","\\3"
    )) %>% plyranges::as_granges()

  return(tbl)

}


#' Fix column names of the SpliceWiz edgeR output.
#'
#' @param x an output of the function [SpliceWiz::ASE_edgeR()].
#'
#' @return a tibble.
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#' // to do
#' }
.fix_splicewiz_res_edger_columns <- function(x){

  x %>% tibble::as_tibble() %>%
    dplyr::rename(c("event_name" = "EventName", "event_type" = "EventType" , "event_region" = "EventRegion", "log2FoldChange" = "logFC", "pvalue" ="PValue", "padj" = "FDR", "delta_psi" = "deltaPSI")) %>%
    dplyr::rename_with(~stringr::str_replace(string = ., pattern = "AvgPSI_(.*)", replacement = "avg_psi_\\1"))
}






#' @title  Categorize diff ASE in to `up` and `down` based on log2FC, p-value and p-adj values.
#' @description Usually fold change 1.5 (log2FC ~0.6) and statistical significance defined by p and/or p.adj values
#' are considered to mark gene (ASE) as differential expressed. However, these cutoffs change based on data quality.
#' This function allows changing these cutoffs and mark genes (ASE) differently expressed.
#' While using this package, most of the time user do not need to call this function explicitly as it is internally called by
#' [parcutils::run_ase_diff_analysis()].
#'
#' @param ase_tibble  a dataframe obtained from \code{SpliceWiz::ASE_*} function.
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
.categorize_diff_ase <- function(ase_tibble,
           log2fc_cutoff =  1,
           pval_cutoff = 0.05,
           padj_cutoff = 0.01,
           add_column_regul = TRUE,
           regul_based_upon = 1
  ){


    # add column 'type' having one of the four values

    # 1. NS -> Not significant
    # 2. p-value -> only pvalue significant but gene is not diff expressed
    # 3. log2FC -> only diff expressed but not pvalue significant
    # 4. p-value&log2FC -> diff expressed and pvalue signficant

    pval_cutoff <- pval_cutoff
    log2fc_cutoff <- log2fc_cutoff

    ase_tibble %<>% dplyr::mutate(signif = dplyr::case_when(dplyr::between(log2FoldChange , -log2fc_cutoff, log2fc_cutoff ) &  pvalue >= pval_cutoff ~ "NS",
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
        ase_tibble %<>%
          dplyr::mutate(regul = dplyr::case_when(log2FoldChange >= log2fc_cutoff & pvalue <= pval_cutoff ~ "Up",
                                                 log2FoldChange <= -1 * log2fc_cutoff & pvalue <= pval_cutoff ~ "Down",
                                                 TRUE ~ "other"))

      }
      # if regul_based_upon == 2

      #   Up -> log2fc >= log2fc_cutoff & padj <= padj_cutoff
      #   Down -> log2fc  <= (-1) * log2fc_cutoff & padj <= padj_cutoff
      #   Other -> remaining genes

      else if(regul_based_upon == 2){
        ase_tibble %<>%
          dplyr::mutate(regul = dplyr::case_when(log2FoldChange >= log2fc_cutoff & padj <= padj_cutoff ~ "Up",
                                                 log2FoldChange <= -1 * log2fc_cutoff & padj <= padj_cutoff ~ "Down",
                                                 TRUE ~ "other"))

      }

      # if regul_based_upon == 3

      #   Up -> log2fc >= log2fc_cutoff & pvalue <= pvalue_cutoff & padj <= padj_cutoff
      #   Down -> log2fc  <= (-1) * log2fc_cutoff pvalue <= pvalue_cutoff & padj <= padj_cutoff
      #   Other -> remaining genes

      else if(regul_based_upon == 3){
        ase_tibble %<>%
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


#' Get Alternate Splice Events based on their differential regulation.
#' @details For a given differential comparison this function returns `up`, `down`, `both`, `other` and `all` ASE.
#'
#' @param x an abject of class "parcutils_ase". This is an output of the function [parcutils::run_ase_diff_analysis()].
#' @param sample_comparisons a character vector denoting  sample comparisons for which ASE to be obtained.
#' @param regulation a character string, default \code{both}. Values can be one of the \code{up}, \code{down}, \code{both}, \code{other}, \code{all}.
#'  + `up` : returns all up regulated ASE.
#'  + `down` : returns all down regulated ASE
#'  + `both` : returns all up and down regulated ASE.
#'  + `other` : returns ASE other than up and down regulated ASE.
#'  + `all` : returns all ASE.
#' @param simplify logical, default FALSE, if TRUE returns result in a dataframe format.
#' @return a list or dataframe.
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#'
#' res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated")
#'
#' # get both up and down regulated genes
#' get_ASE_by_regulation(x = res, sample_comparisons = c("A_VS_B")) %>% head()
#'
#' # get up genes only
#' get_ASE_by_regulation(x = res, sample_comparisons = c("A_VS_B") , regul = "up") %>% head()
#'
#' # get down genes only
#' get_ASE_by_regulation(x = res, sample_comparisons = c("A_VS_B") , regul = "down") %>% head()
#'
#' # get genes other than up and down
#' get_ASE_by_regulation(x = res, sample_comparisons = c("A_VS_B") , regul = "other") %>% head()
#'
#' # Simplify output for multiple sample comparisons
#' get_ASE_by_regulation(x = res, sample_comparisons = res$de_comparisons, simplify = TRUE, regul= "up")
#'
#'
#' # get genesets by regulation. It uses sample comparison and regulation to name each output geneset.
#'
#' get_ASEsets_by_regulation(x = res, sample_comparisons = "A_VS_B", regul = "both")
#'
get_ASE_by_regulation <- function(x, sample_comparisons , regulation = "both" , simplify = FALSE  ) {

  # validate x.
  stopifnot("x must be an object of class 'parcutils_ase'. Usually x is derived by parcutils::run_deseq_analysis()." = is(x, "parcutils_ase"))

  # validate sample comparison.
  stopifnot("sample_comparisons must be a character vector" = is.character(sample_comparisons) & length(sample_comparisons) >= 1)

  # validate regulation.
  match.arg(regulation , choices = c("up","down" , "both", "other" ,"all"), several.ok = FALSE)

  # validate simplify
  stopifnot("simplify must be a logical value" = is.logical(simplify))

  # if length of  sample_comparisons > 1, use this function recursively
  if(length(sample_comparisons) > 1){
    rslt <- purrr::map(sample_comparisons,
                       ~get_ASE_by_regulation(x = x,sample_comparisons = ..1,regulation = regulation , simplify = simplify))
    names(rslt) <- sample_comparisons
    return(rslt)
  }

  # split by column regul.
  # NOTE: do not group by column name instead use index. The reason is because if the column name of column 'regul' change in future it will break this code.

  ase_by_comp <- x$res_ase_diff_tibble_annot[[sample_comparisons]] %>%
    dplyr::select(1, dplyr::last_col())

  ase_by_comp <- ase_by_comp %>%
    dplyr::group_by_at(2)

  grp_keys <- ase_by_comp %>%
    dplyr::group_keys() %>%
    dplyr::pull(1)

  ase_by_comp <- ase_by_comp %>%
    dplyr::group_split()
  names(ase_by_comp) <- grp_keys

  # Convert names and regulation to lower case before comparison
  names(ase_by_comp)  <- tolower(names(ase_by_comp))
  regulation = tolower(regulation)


  if(regulation == "both"){
    # select up and/or down
    rslt <- ase_by_comp[names(ase_by_comp) %in% c("up","down")]  %>% as.list() %>% dplyr::bind_rows()
  } else if(regulation == "all"){
    # select all
    rslt <- ase_by_comp  %>% as.list() %>% dplyr::bind_rows()
  } else{
    # select other
    rslt <- ase_by_comp[names(ase_by_comp) %in% regulation]  %>% as.list() %>% dplyr::bind_rows()
  }

  # convert factor
  if(nrow(rslt) == 0){
    rslt <- factor(levels = c("up","down","other"))
  } else{
    rslt <- factor(rslt[[2]] %>% tolower(), levels = c("up","down","other"))  %>%
      purrr::set_names(rslt[[1]])
  }

  if(simplify){
    rslt %<>% tibble::enframe() %>%
      dplyr::rename(c("regul" = "value" , "ase" = "name")) %>%
      dplyr::mutate(sample_comparisons = sample_comparisons) %>%
      dplyr::select(dplyr::all_of(c("sample_comparisons", "ase", "regul")))
  }
  return(rslt)

}



#' @rdname get_ASE_by_regulation
#' @export
get_ASEsets_by_regulation <- function(x, sample_comparisons, regulation = "both" ){


  ase <- get_ASE_by_regulation(x, sample_comparisons = sample_comparisons,
                                   regulation = regulation,
                                   simplify = TRUE)

  ase_sets <- dplyr::bind_rows(ase) %>%
    dplyr::mutate(ase_set_name = stringr::str_c(sample_comparisons, regul , sep = "_")) %>%
    named_group_split(ase_set_name) %>%
    purrr::map(~ ..1 %>% dplyr::pull(ase))

  return(ase_sets)
}


#' Generate a volcano plot for diff. ASE.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_ase_diff_analysis()].
#' @param sample_comparison a character string denoting a valid differential ASE comparison. Possible comparisons can be found from x$de_comparisons.
#' @param log2fc_cutoff a numeric value denoting logFC cutoff for volcano plot. Default 1.
#' @param pval_cutoff a numeric value denoting pvalue cutoff for volcano plot. Default 0.05.
#' @param ASE_to_display a character vector of the ASE to display in volcano plot, default NULL, displays non overlapping ASE.
#' @param lab_size a numeric value, default 3, denoting size of the labels.
#' @param point_size a numeric value, default 1, denoting size of the points.
#' @param repair_ASE logical, default TRUE, indicating whether to repair event names. See details.
#' @param col_up a character string, default \code{#a40000}, denoting valid color code for up regulated ASE.
#' @param col_down a character string, default \code{#007e2f}, denoting valid color code for up regulated ASE.
#' @param col_other a character string, default "grey", denoting valid color code for other than up and down regulated ASE.
#' @param facet_event_type logical, default TRUE, denoting whether to facet plot by event type.
#' @param facte_nrow pass to the argument `nrow` of [ggplot2::facet_wrap()].
#' @param facet_ncol pass to the argument `ncol` of [ggplot2::facet_wrap()].
#' @param facet_scale pass to the argument `scale` of [ggplot2::facet_wrap()].
#' @param facet_shrink pass to the argument `shrink` of [ggplot2::facet_wrap()].
#' @param facet_labeller pass to the argument `labeller` of [ggplot2::facet_wrap()].
#' @param facet_as.table pass to the argument `as.table` of [ggplot2::facet_wrap()].
#' @param facet_drop pass to the argument `drop` of [ggplot2::facet_wrap()].
#' @param facet_dir pass to the argument `dir` of [ggplot2::facet_wrap()].
#' @param facet_strip.position pass to the argument `position` of [ggplot2::facet_wrap()].
#' @param ... other parameters to be passed to [EnhancedVolcano::EnhancedVolcano()].
#' @details
#' + repair_ASE: Internally event names are taken from SpliceWiz. When repair_genes is set to TRUE the string corresponding to gene symbol will be extracted.
#'  if one gene symbol assinged to more than one events, even names will be used instead of gene symbol (usually the case). This is useful when gene names to be revealed in the volcano plot.
#' @return a ggplot2.
#' @export
#'
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' res <- run_ase_diff_analysis(x = se, test_factor = "treatment", test_nom = "A" ,test_denom = "B",  IRmode ="annotated",  cutoff_lfc = 0.6, cutoff_pval = 1, regul_based_upon = 1)
#' get_ase_volcano_plot(res,  sample_comparison = res$de_comparisons[1], pval_cutoff = 1) %>% print()
#'
#' # repair_ASE = FALSE
#' get_ase_volcano_plot(res,  sample_comparison = res$de_comparisons[1], pval_cutoff = 1,repair_ASE = FALSE) %>% print()
#'
#' # ASE_to_display = "", don't show any label.
#' get_ase_volcano_plot(res,  sample_comparison = res$de_comparisons[1], pval_cutoff = 1,ASE_to_display = "") %>% print()
get_ase_volcano_plot <- function(x,
         sample_comparison,
         log2fc_cutoff = 1,
         pval_cutoff = 0.05,
         ASE_to_display = NULL,
         lab_size = 3,
         point_size = 1,
         repair_ASE = TRUE,
         col_up = "#a40000",
         col_down =  "#007e2f",
         col_other = "grey",
         facet_event_type = TRUE,
         facte_nrow = NULL,
         facet_ncol = NULL,
         facet_scale = "fixed",
         facet_shrink = TRUE,
         facet_labeller = "label_value",
         facet_as.table = TRUE,
         facet_drop = TRUE,
         facet_dir = "h",
         facet_strip.position = "top",
         ...){

  # validate x
  .validate_parcutils_ase_obj(x)

  # validate sample_comparison

  stopifnot("sample_comparison must be a character string" = is.character(sample_comparison) & length(sample_comparison) == 1)

  # validate genes_to_display

  stopifnot("ASE_to_display must be a character vector" = is.character(ASE_to_display) | is.null(ASE_to_display))


  # prepare volcano plots

  # filter by sample_comparison
  volcano_data <- x$res_ase_diff_tibble[[sample_comparison]] %>%dplyr::select(event_name, log2FoldChange , pvalue, event_type)


  # Fix ASE names. It gives gene names if it's unique. Otherwise event names will be returned.

  if(repair_ASE){
    an <- volcano_data$event_name
    gsym <- get_genes_from_event_name(an)
    an_new <- dplyr::if_else(duplicated(gsym) , an , gsym)
    volcano_data$event_name <- an_new
  }

  # generate plot

  pp <- parcutils::EnhancedVolcano2(toptable = volcano_data,
                                    x = "log2FoldChange" ,
                                    y = "pvalue",
                                    lab = volcano_data$event_name,
                                    selectLab = ASE_to_display,
                                    FCcutoff = log2fc_cutoff,
                                    pCutoff = pval_cutoff,
                                    labSize = lab_size,
                                    pointSize =point_size,
                                    col_by_regul = TRUE,
                                    col_up = col_up,
                                    col_down = col_down,
                                    col_other = col_other,
                                    title = sample_comparison)#  , ...)


  if(facet_event_type){
    pp <- pp + ggplot2::facet_wrap(~event_type,
                             nrow = facte_nrow,
                             scales = facet_scale, ncol = facet_ncol, shrink = facet_shrink,labeller = facet_labeller,as.table = facet_as.table,
                             drop = facet_drop, dir = facet_dir,strip.position =  facet_strip.position)
  }

  return(pp)


}

# validate parcutils_ase object
.validate_parcutils_ase_obj <- function(x){
    stopifnot("x must be an object of class 'parcutils_ase'. Usually x is derived by parcutils::run_ase_diff_analysis()." = is(x, "parcutils_ase"))
}
