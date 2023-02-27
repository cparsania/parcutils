#' Generate a barplot showing counts of splice events
#'
#' @param x an object of class \code{NxtSE}
#'
#' @return ggplot.
#' @export
#'
#' @examples
#' \dontrun{
#' // to do
#' }
#'
#' @keywords internal
spliceWiz_plot_event_counts <- function(x){

  stopifnot("x must be the class of NxtSE" = is(x, "NxtSE"))
  gp <- se.filtered@elementMetadata %>%
    tibble::as_tibble() %>%
    dplyr::select(EventName, EventType) %>%
    dplyr::group_by(EventType) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(EventType = forcats::fct_reorder(EventType, -n)) %>%
    ggplot(aes(x = EventType, y = n)) +
    geom_bar(position ="dodge", stat = 'identity') +
    geom_text(aes(label=n), inherit.aes = TRUE, position = position_dodge(width = 0.9),size = 4, vjust=-0.25) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 15)) +
    ggplot2::ylab("Counts") +
    ggplot2::xlab("EventType") +
    ggeasy::easy_rotate_labels(which = "x") +
    theme(text = element_text(size = 15))
  return(gp)

}



#' Get diff event names from SpliceWiZ edger diff results
#'
#' @param res_edger an output of [SpliceWiz::ASE_edgeR()]
#' @param event_type a character string from one of the followings: SE,A5SS, A3SS, MXE, I,AFE, ALE. Default SE
#' @param padj_cutoff a double denoting cutoff value of padj (FDR). Default 0.05
#' @param log2fc_cutoff a double denoting cutoff value of log2FC (logFC). Default 0.6
#'
#' @return a character vector
#' @export
#'
#' @examples
#'
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' require("edgeR")
#' res_edgeR <- SpliceWiz::ASE_edgeR(se, "treatment", "A", "B", useQL = FALSE)
#' res_edgeR_QL <- SpliceWiz::ASE_edgeR(se, "treatment", "A", "B", useQL = TRUE)
#' get_edger_diff_ase_event_names(res_edgeR,padj_cutoff= 1)

get_edger_diff_ase_event_names <- function(res_edger,
                                           event_type = "SE",
                                           padj_cutoff= 0.05,
                                           log2fc_cutoff = 0.6){

  res_edger %>%
    dplyr::filter(EventType == !!event_type) %>%
    dplyr::filter((logFC <= -(log2fc_cutoff) | logFC >= log2fc_cutoff) &
                    FDR <= padj_cutoff) %>% dplyr::pull(EventName)


}

#' Get an input matrix for diff ASE heatmap
#'
#' @param se an object of class NxtSE.
#' @param event_names a character vector denoting valid ASE names to plot in heatmap.
#' @param samples a character vector denoting valid sample names to plot in the heatmap.
#' @param replicates logical, default TRUE, denoting whether to show sample replicates or not.
#' @param method a character string denoting a method for data value. Can be one of the "PSI", "logit" or "Z-Score"
#' @param column_condtion a name of column in se storing condition. Deafault: "condition"
#'
#' @return a data matrix
#' @export
#' @keywords internal
#' @examples
#' se <- SpliceWiz::SpliceWiz_example_NxtSE(novelSplicing = TRUE)
#' SpliceWiz::colData(se)$treatment <- rep(c("A", "B"), each = 3)
#' SpliceWiz::colData(se)$replicate <- rep(c("P","Q","R"), 2)
#' res_edgeR <- SpliceWiz::ASE_edgeR(se, "treatment", "A", "B", useQL = FALSE)
#' event_names = get_edger_diff_ase_event_names(res_edgeR,   padj_cutoff = 1)
#' .get_aes_heatmap_data(se, event_names , samples = c("A", "B"), column_condition ="treatment")
.get_aes_heatmap_data <- function(se , event_names,samples,replicates=TRUE, method ="PSI", column_condition ="condition"){

  samples_repli <- SpliceWiz::colData(se) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "samples") %>%
    dplyr::filter(!!rlang::sym(column_condition) %in% !!samples) %>%
    dplyr::pull(samples)

  data_matrix <- SpliceWiz::makeMatrix(se = se,
                                       event_list = event_names,
                                       sample_list = samples_repli,
                                       method = method)

  return(data_matrix)
}



#' Generate a heatmap for SpliceWiz splice events
#'
#' @param se an object of class NxtSE.
#' @param event_names a character vector denoting valid spliceWiz event names
#' @param samples a character vector denoting valid sample names
#' @param column_condition a name of column in se storing condition. Deafault: "condition"
#' @param method one of the character strings: "Z-score", "PSI" ,"logit"
#' @param show_replicates logical, default TRUE, denoting whether to show sample replicates or not.
#' @param show_row_names logical, whether to show row names in the heatmap.
#' @param cluster_rows logical, whether to clusters rown in the heatmap
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
#' res_edgeR <- SpliceWiz::ASE_edgeR(se, "treatment", "A", "B", useQL = FALSE)
#' event_names = get_edger_diff_ase_event_names(res_edgeR,   padj_cutoff = 1)
#' get_diff_ase_heatmap(se, event_names = event_names, samples = c("A" ,"B"), column_condition    = "treatment")
#' get_diff_ase_heatmap(se, event_names = event_names, samples = c("A" ,"B"), column_condition    = "treatment",method = "Z-score")
get_diff_ase_heatmap <- function(se,
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


  dat <- .get_aes_heatmap_data(se = se,
                               event_names = event_names ,
                               column_condition = column_condition,
                               samples = samples,method = method)

  dat <- dat %>% as.data.frame() %>% tibble::rownames_to_column() %>%
    tibble::as_tibble() %>%
    TidyWrappers::tbl_remove_rows_NA_any() %>% as.data.frame() %>%
    tibble::column_to_rownames("rowname")

  ComplexHeatmap::Heatmap(dat, show_row_names = show_row_names,
                          cluster_rows = cluster_rows,
                          show_column_dend = show_column_dend,
                          cluster_columns = cluster_columns,
                          show_row_dend = show_row_dend, ...)


}



#' Get gene names from spliceWiz event names
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





