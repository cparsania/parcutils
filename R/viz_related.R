
#' Sort order of box obtained from geom_boxplot
#'
#' @param x ggplot2 with atleast one layer of geom_boxplot
#' @param decreasing logical, default T, indicating whether the box are arraged in increasing or decreasing order.
#'
#' @return a ggplot.
#' @export
#' @importFrom purrr map_lgl
#' @importFrom rlang expr
#' @importFrom dplyr summarise
#' @importFrom dplyr arrange
#' @importFrom dplyr pull
#' @importFrom dplyr group_by
#' @importFrom dplyr mutate
#' @examples
#' library(ggplot2)
#' dd <- iris %>%
#'   tidyr::pivot_longer(cols = Sepal.Length:Petal.Width)
#'
#' # single layer
#' ss <- dd %>%
#'   ggplot2::ggplot() +
#'   geom_boxplot(aes(x = name, y = value))
#' sort_geom_box(x = ss,decreasing = T)
#' sort_geom_box(x = ss,decreasing = F)
#'
#' # multiple layer
#' ss <- dd %>%
#'   ggplot2::ggplot() +
#'   geom_point(aes(y = value , x = name)) +
#'   geom_boxplot(aes(x = name, y = value))
#' sort_geom_box(x = ss,decreasing = T)
#' sort_geom_box(x = ss,decreasing = F)
#'
sort_geom_box <- function(x, decreasing = T){
  #x <- ss
  plt <- x

  stopifnot(inherits(plt , "ggplot"))

  # get data from ggplot object
  dat <- plt$data

  #which layer is "GeomBoxplot"
  is_ggbox <- purrr::map_lgl(plt$layers, ~inherits(..1$geom,"GeomBoxplot"))

  if(all(!is_ggbox)){
    stop("'x'  must contain atleast one layer of 'geom_boxplot'")
  }

  is_ggbox_index <- which(is_ggbox == T)

  oo <- plt$layers[[is_ggbox_index]]

  # get var names of x and y
  x <- oo$mapping$x
  y <- oo$mapping$y

  # order -> decrease or increase
  if(decreasing){
    dec <- rlang::expr(-median)
  }else {
    dec <- rlang::expr(median)
  }


  # get levels of x sorted by median
  new_lvls <- dat %>%
    dplyr::group_by(!!x) %>%
    dplyr::summarise(median =median(!!y)) %>%
    dplyr::arrange(!!dec) %>%
    dplyr::pull(!!x) %>%
    as.character()

  # set new levels in the data
  dat <- dat %>%
    dplyr::mutate(!!x :=  forcats::fct_relevel(!!x, new_lvls))

  # update data into plot
  plt$data <- dat

  return(plt)
}


# volcano plot

#' Wrapper around [EnhancedVolcano::EnhancedVolcano()]
#'
#' @param toptable \code{toptable} [EnhancedVolcano::EnhancedVolcano()]
#' @param lab \code{lab} [EnhancedVolcano::EnhancedVolcano()]
#' @param x \code{x} [EnhancedVolcano::EnhancedVolcano()]
#' @param y \code{y} [EnhancedVolcano::EnhancedVolcano()]
#' @param pCutoff \code{pCutoff} [EnhancedVolcano::EnhancedVolcano()]
#' @param FCcutoff \code{FCcutoff} [EnhancedVolcano::EnhancedVolcano()]
#' @param col_by_regul lgl, default T, whether to color variable by gene regulation - Up, Down, Other
#' @param ... Other parameters pass to [EnhancedVolcano::EnhancedVolcano()]
#' @param col_up a character string denoting color for up genes, works only when col_by_regul = T
#' @param col_down  a character string denoting color for down genes, works only when col_by_regul = T.
#' @param col_others a character string denoting color for other genes, works only when col_by_regul = T.
#'
#' @return a volcano plot
#' @export
EnhancedVolcano2 <- function(toptable, lab, x, y, pCutoff = 10e-4,FCcutoff = 1.5,col_by_regul = T,
                             col_up = "#b2182b" ,col_down = "#2166ac", col_others = "#e0e0e0" ,...){

  if(col_by_regul) {

    keyvals <- ifelse(
      toptable[[x]] <= -FCcutoff  & toptable[[y]] <= pCutoff , col_down,
      ifelse(toptable[[x]] >= FCcutoff &  toptable[[y]] <= pCutoff, col_up,
             col_others))
    keyvals[is.na(keyvals)] <- col_others
    names(keyvals)[keyvals == col_down ] <- 'Down'
    names(keyvals)[keyvals == col_others] <- 'Other'
    names(keyvals)[keyvals == col_up] <- 'Up'

    plot  <- EnhancedVolcano::EnhancedVolcano(toptable = toptable, lab = lab,x = x,y = y,
                                              pCutoff = pCutoff, FCcutoff = FCcutoff, colCustom  = keyvals, ...)
  } else {
    plot  <- EnhancedVolcano::EnhancedVolcano(toptable = toptable, lab = lab,x = x,y = y,
                                              pCutoff = pCutoff, FCcutoff = FCcutoff, ...)
  }
  plot
}







# visualize relative position of regions related to reference region



#' Generate segment plot to visualize genomic regions related to reference regions.
#'
#' @param query an object of [GenomicRanges::granges()] to be visualized relative to reference
#' @param reference an object of [GenomicRanges::granges()] against which \code{query} regions should be mapped.
#'
#' @return a list containing segment plot and plot data
#' @export
#'
#' @examples
#' \dontrun{
#' ref_file <- system.file("extdata", "Homo_sapiens.GRCh38.101.gtf.gz" ,package = "parcutils")
#' ref_data <- rtracklayer::import(ref_file)  %>%
#'  tibble::as_tibble() %>%
#'  dplyr::mutate(seqnames = stringr::str_c("chr",seqnames ,sep = "")) %>% plyranges::as_granges()
#'  ref_data_gr <-  ref_data %>%
#' dplyr::filter(type == "gene" & gene_biotype == "protein_coding")
#' query_regions_gr <- ref_data %>% dplyr::filter(type == "Selenocysteine")  %>% sample(5000)
#'
#' oo <- plot_regions_relative_to_reference(query = query_regions_gr, reference  =  ref_data_gr)
#'  oo$plot
#' oo$data %>% ggplot2::ggplot(aes(x = q_relt_start)) + geom_density()
#'}
plot_regions_relative_to_reference <- function(query , reference){

  # query <- query_regions_gr
  # reference  <-  ref_data_gr


  # check object type. Both query and reference must be of the class gr

  if(!any("GRanges" %in% attr(query, "class") )){
    stop("'query' must be of an object of class 'GRanges'")
  }

  if(!any("GRanges" %in% attr(reference, "class") )){
    stop("'reference' must be of an object of class 'GRanges'")
  }

  ##
  # remove duplicates
  ##

  query <- query %>% unique()
  reference <- reference %>% unique()

  ##
  # add unique identifier
  ##

  # id for query
  query_id_col_name <- "query_id"
  query <- query %>%
    tibble::as_tibble() %>%
    dplyr::mutate(!!query_id_col_name := stringr::str_c(seqnames , start, end , strand , sep = "|")) %>%
    plyranges::as_granges()

  # id for reference

  ref_id_col_name <- "reference_id"
  reference <- reference %>%
    tibble::as_tibble() %>%
    dplyr::mutate(!!ref_id_col_name := stringr::str_c(seqnames , start, end , strand , sep = "|")) %>%
    plyranges::as_granges()

  ##
  # find reference for each query region.
  ##

  # query is within reference

  if(TRUE){
    q_r_combined <- query %>%
      plyranges::join_overlap_left_within_directed(reference)
  }

  # mapping statistics

  all_queries <- query$query_id
  all_reference <- reference$reference_id

  mapped_id_table <- q_r_combined %>%
    tibble::as_tibble() %>%
    dplyr::select(!!query_id_col_name , !!ref_id_col_name) %>%
    # if no reference mapped to query then it will be NA.
    tidyr::drop_na() %>%
    unique()

  mapping_stats <-  tibble::tibble(all_queries = all_queries) %>%
    dplyr::mutate(times_mapped = purrr::map_dbl(all_queries , ~ which(..1 ==  mapped_id_table$query_id) %>% length()))

  # prepare final table for plot

  for_plot <- mapped_id_table %>%
    tidyr::separate(col = query_id , into = c("qseqname","qstart", "qend", "qstrand"), sep = "\\|", convert = T ,remove = F) %>%
    tidyr::separate(col = reference_id , into = c("rseqname","rstart", "rend", "rstrand") , sep = "\\|", convert = T, remove = F) %>%

    # add width
    dplyr::mutate(qwidth = qend - qstart + 1 , rwidth = rend - rstart + 1) %>%

    # add absolute start and end
    dplyr::mutate(q_abs_start = 1 + abs(qstart - rstart), q_abs_end = q_abs_start + (qend - qstart)) %>%

    # add relative start and end. Strand must be considered
    dplyr::mutate(q_relt_start = dplyr::if_else(qstrand == "-", 1- (q_abs_start/rwidth)  , q_abs_start/rwidth) ,
                  q_relt_end =  dplyr::if_else(qstrand == "-" , 1 - (q_abs_end/rwidth) ,  q_abs_end/rwidth )) %>%

    # flip start and end after strand considerd. This will make sure that start <  end
    dplyr::mutate(q_relt_start_temp = dplyr::if_else(qstrand == "-",  q_relt_end , q_relt_start) ,
                  q_relt_end_temp = dplyr::if_else(qstrand == "-",  q_relt_start , q_relt_end)) %>%

    # remove temp cols
    dplyr::mutate(q_relt_start = q_relt_start_temp , q_relt_end = q_relt_end_temp) %>%
    dplyr::select(-q_relt_start_temp,-q_relt_end_temp)


  # add density
  if(TRUE) {
    dens <- density(for_plot$q_relt_start)

    #approx density (y) for given (x)
    af  <-  approxfun(dens$x , dens$y)

    density_col_name <- "dens"
    for_plot %<>% dplyr::mutate(!!density_col_name := af(q_relt_start))

  }

  # sort by start postion of query region
  if(TRUE){
    for_plot %<>%
      dplyr::mutate(id = paste0(dplyr::row_number())) %>%
      dplyr::mutate(id = forcats::fct_reorder(id, q_relt_start))
  }

  # plot
  gp <- for_plot %>%
    ggplot2::ggplot(ggplot2::aes(x = q_relt_start, y = id)) +
    ggplot2::geom_segment(ggplot2::aes(xend = q_relt_end, yend = id, col = dens)) +
    ggplot2::scale_color_viridis_c() +
    xlim(0, 1) +
    theme_bw() +
    labs(color = "Density") +
    theme(text = element_text(size = 15),
          axis.text = element_text(size = 15), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), legend.text = element_text(size = 12)) +
    xlab("Relative query start") +
    ylab("Query regions")


  rtrn <- list(plot= gp, data = for_plot)

}



#' Generate a PCA plot.
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting samples to plot in PCA plot, default \code{NULL}. If set to NULL all samples are accounted.
#' @param genes a character vector denoting genes to consider for PCA plot, default \code{NULL}. If set to NULL all genes are accounted.
#' @param circle_size a numeric value,  default  10, denoting size of the circles in PCA plot.
#' @param show_replicates a logical, degfault FALSE, denoting whether to label replicates in the PCA plot.
#'
#' @return an object of ggplot2.
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
#'                          group_denominator = c("control"),
#'                          column_samples = c("control_rep1", "treat1_rep1", "treat2_rep1", "control_rep2", "treat1_rep2", "treat2_rep2", "control_rep3", "treat1_rep3", "treat2_rep3"))
#'
#'
#'  get_pca_plot(x = res) %>% print()
#'
#' # label replicates
#'
#' get_pca_plot(x = res, label_replicates =  T)
#'
get_pca_plot <- function(x, samples = NULL, genes = NULL, circle_size = 10, label_replicates = FALSE){

  # validate x
  parcutils:::validata_parcutils_obj(x)

  # validate samples
  stopifnot("'samples' must be a character vector or NULL." = is.character(samples) | is.null(samples))

  # validate genes
  stopifnot("'genes' must be a character vector or NULL." = is.character(genes) | is.null(genes))

  # validate show_replicates
  stopifnot("'label_replicates' must be a logical." = is.logical(label_replicates))

  # get gene expression matrix

  norm_expr_mat <- parcutils::get_normalised_expression_matrix(x = x, samples = samples, genes = genes , summarise_replicates = FALSE)

  # convert log2

  norm_expr_mat_log2 <- norm_expr_mat %>% dplyr::mutate_if(is.numeric, ~log2(. + 0.1))

  # get column name for gene id

  gene_id_col = norm_expr_mat_log2 %>% colnames() %>% .[1]


  # prepare sample for PCA.

  pca_input <- norm_expr_mat_log2  %>%
    tidyr::pivot_longer(names_to = "sample" ,values_to  = "value" , -!!rlang::sym(gene_id_col)) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(gene_id_col), values_from = !!rlang::sym("value")) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("sample")

  pr_comps <- prcomp(pca_input,scale = T)

  pr_comps_tbl <- pr_comps %>%
    .$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("samples") %>%
    tibble::as_tibble()

  # get prop of variance

  pc_prop_of_var <- summary(pr_comps)$importance %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "feature_type") %>%
    tidyr::gather("pc" , "value" , -feature_type) %>%
    dplyr::filter(feature_type == "Proportion of Variance") %>%
    dplyr::select(-dplyr::all_of("feature_type")) %>%
    dplyr::mutate(value = round(value * 100 ,1)) %>% ## convert to percentage
    dplyr::mutate(with_var = paste(pc ," (" , value ,"%",")" , sep = "")) %>%
    dplyr::select(-dplyr::all_of("value"))

  # return named vector

  pc_prop_of_var <- pc_prop_of_var %>%
    dplyr::pull(with_var) %>%
    rlang::set_names( pc_prop_of_var$pc)

  # prepare sample_info from x

  sample_info <- group_replicates_by_sample(x = x)

  # add sample info

  pr_comps_tbl %<>% dplyr::left_join(sample_info , by = c("samples"))

  # select which PCs to plot

  x_pc = rlang::quo(PC1)
  y_pc = rlang::quo(PC2)

  # pca plot  ----

  pca_plot <- pr_comps_tbl %>%
    ggplot2::ggplot(ggplot2::aes(x = !!x_pc, y = !!y_pc)) +
    ggplot2::geom_point(ggplot2::aes(fill = groups),size = circle_size, pch=21)

  ## show replicates
  if(label_replicates){
    pca_plot <- pca_plot + ggrepel::geom_label_repel(ggplot2::aes(label = samples),max.overlaps = 100000)
  }

  # alter theme

  pca_plot <- pca_plot + ggplot2::theme_bw()+
    ggplot2::theme(text = ggplot2::element_text(size = 15,colour = "black"))+
    ggplot2::xlab(pc_prop_of_var[[rlang::as_label(x_pc)]]) +
    ggplot2::ylab(pc_prop_of_var[[rlang::as_label(y_pc)]])

  # pca_plot

  pca_plot <- pca_plot +
    ggplot2::scale_fill_manual(values = MetBrewer::met.brewer(n = 6 ,
                                                              name = "Austria")) +
    ggplot2::guides(fill = ggplot2::guide_legend("Samples"))



  return(pca_plot)

}


#' Generate a box plot
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting samples to plot in boxplot, default \code{NULL}. If set to NULL all samples are accounted.
#' @param genes a character vector denoting genes to consider for boxplot, default \code{NULL}. If set to NULL all genes are accounted.
#' @param group_replicates logical, default \code{FALSE}, whether to group replicates in individual plots.
#' @param convert_log2 logical, default \code{TRUE}, whether to plot log2(noralised expression + 1) or not.
#'
#' @return an object of class ggplot2
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
#'                          group_denominator = c("control"),
#'                          column_samples = c("control_rep1", "treat1_rep1", "treat2_rep1", "control_rep2", "treat1_rep2", "treat2_rep2", "control_rep3", "treat1_rep3", "treat2_rep3"))
#'
#' get_gene_expression_box_plot(x = res ) %>% print()
#'
#' get_gene_expression_box_plot(x = res , group_replicates = T ) %>% print()
#'
get_gene_expression_box_plot <- function(x, samples = NULL, genes = NULL,
                                         group_replicates = F,
                                         convert_log2 = T){


  # validate x
  parcutils:::validata_parcutils_obj(x)

  # validate samples
  stopifnot("'samples' must be a character vector or NULL." = is.character(samples) | is.null(samples))

  # validate genes
  stopifnot("'genes' must be a character vector or NULL." = is.character(genes) | is.null(genes))

  # validate show_replicates
  stopifnot("'group_replicates' must be a logical." = is.logical(group_replicates))


  # get norm counts
  norm_counts <- parcutils::get_normalised_expression_matrix(x = x,
                                                             samples = samples,
                                                             genes = genes,
                                                             summarise_replicates = FALSE)

  column_gene_id <-  colnames(norm_counts)[1]

  # convert log2

  if(convert_log2){
    norm_counts <- norm_counts %>%
      dplyr::mutate_if(is.numeric , ~log2(. + 1))
  }

  # long format
  norm_counts_long <- norm_counts %>%
    tidyr::pivot_longer(cols = -!!rlang::sym(column_gene_id),
                        names_to = "samples",
                        values_to ="value")

  # group replicates by samples

  rep_grps <- parcutils::group_replicates_by_sample(x)

  # add groups

  norm_counts_long %<>% dplyr::left_join(rep_grps , by = ("samples"))

  if(!is.null(samples)){

    # order samples in replicate groups by in order of user input samples.

    norm_counts_long$groups <- forcats::fct_relevel(norm_counts_long$groups , samples)
    norm_counts_long$samples <- factor(norm_counts_long$samples ,
                                       levels = norm_counts_long$samples %>% unique())
  }


  # generate box plot

  box_plots <- norm_counts_long %>%
    ggplot2::ggplot() +
    ggplot2::geom_boxplot(ggplot2::aes(x = samples, y = value, fill =  groups), alpha = 0.9)

  ## theme

  box_plots  <- box_plots + ggplot2::theme_bw()+
    ggplot2::theme(text = ggplot2::element_text(size = 15) ,
                   axis.text.x  = ggplot2::element_text(angle = 90))


  ## col

  box_plots <- box_plots +
    ggplot2::scale_fill_manual(values = MetBrewer::met.brewer(n = 6 , name = "Austria"))  +
    ggplot2::guides(fill = ggplot2::guide_legend("Samples"))


  ## axis lab

  if(convert_log2){
    box_plots <- box_plots + ggplot2::ylab("Log2(Normalised count + 1)")
  } else{
    box_plots <- box_plots + ggplot2::ylab("Normalised counts")
  }

  ## group replicates
  if(group_replicates){
    box_plots <- box_plots + ggplot2::facet_wrap(~groups , scales = "free_x")
  }

  return(box_plots)
}



#' Get a scatter plot showing correlation between replicates.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param samples a character vector denoting samples to plot in scatter plot, default \code{NULL}. If set to NULL all samples are accounted.
#' @param genes a character vector denoting genes to consider in scatter plot, default \code{NULL}. If set to NULL all genes are accounted.#'
#' @return a named list.
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
#'                          group_denominator = c("control"),
#'                          column_samples = c("control_rep1", "treat1_rep1", "treat2_rep1", "control_rep2", "treat1_rep2", "treat2_rep2", "control_rep3", "treat1_rep3", "treat2_rep3"))
#'
#' cp <- get_pairwise_corr_plot(res)
#'
#'
#' names(cp) %>% print()
#'
#' cp[1] %>% print
#'
#' cp[2] %>% print
#'
#' cp[3] %>% print
#'
get_pairwise_corr_plot  <- function(x , samples = NULL, genes = NULL){

  # validate x
  parcutils:::validata_parcutils_obj(x)

  # validate samples
  stopifnot("'samples' must be a character vector or NULL." = is.character(samples) | is.null(samples))

  # validate genes
  stopifnot("'genes' must be a character vector or NULL." = is.character(genes) | is.null(genes))

  # get expression values
  norm_expr_mat <- parcutils::get_normalised_expression_matrix(x = x, samples = samples, genes = genes, summarise_replicates = F)

  # convert log2
  if(TRUE){
    norm_expr_mat <- norm_expr_mat %>% dplyr::mutate_if(is.numeric, ~ log2(. + 1))
  }

  # replicate groups
  rep_grps <- get_replicates_by_sample_list(x)

  # filter by user supplied samples
  if(!is.null(samples)){
    rep_grps <- rep_grps[samples]
  }

  # plot
  replicate_corr_plts <- purrr::map(rep_grps , ~ GGally::ggpairs(norm_expr_mat,
                                                                 columns = ..1) +
                                      ggplot2::theme_bw() +
                                      ggplot2::theme(text = ggplot2::element_text(size = 20)))

  return(replicate_corr_plts)

}


#' Group replicates by sample - returns list
#' @rdname get_replicates_by_sample
#'
#' @keywords internal
#' @export
#' @examples
#' \dontrun{
#' // TO DO
#' }
get_replicates_by_sample_list <- function(x){

  # validate x
  parcutils:::validata_parcutils_obj(x)

  y <- parcutils::group_replicates_by_sample(x)

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




#' Generate a volcano plot.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#' @param sample_comparison a character string denoting a valid differatnial gene comparison. Possible comparisons can be found from x$comp.
#' @param log2fc_cutoff a numeric value, default 1.
#' @param pval_cutoff a numeric value, default 0.05.
#' @param genes_to_display a character vector of the genes to display in volcano plot, default NULL, displays non overlapping genes.
#' @param lab_size a numeric value, default 3, denoting size of the lables.
#' @param point_size a numeric value, default 1, denoting size of the points
#' @param col_up a character string, default "a40000", denoting valid color code for up regulated genes.
#' @param col_down a character string, default "007e2f", denoting valid color code for down regulated genes.
#' @param col_other a character string, default "grey", denoting valid color code for other than up and down regulated genes.
#' @param ... other parameters to be passed to [EnhancedVolcano::EnhancedVolcano()]
#'
#' @return ggplot
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
#'                          group_denominator = c("control"),
#'                          column_samples = c("control_rep1", "treat1_rep1", "treat2_rep1", "control_rep2", "treat1_rep2", "treat2_rep2", "control_rep3", "treat1_rep3", "treat2_rep3"))
#'
#'
#' get_volcano_plot(res,  sample_comparison = res$comp[1]) %>% print()
#'
#' get_volcano_plot(res,  sample_comparison = res$comp[2]) %>% print()
#'
get_volcano_plot <- function(x,
                             sample_comparison,
                             log2fc_cutoff = 1,
                             pval_cutoff = 0.05,
                             genes_to_display = NULL,
                             lab_size = 3,
                             point_size = 1,
                             col_up = "#a40000",
                             col_down =  "#007e2f",
                             col_other = "grey",...){

  # validate x
  parcutils:::validata_parcutils_obj(x)

  # validate sample_comparison

  stopifnot("sample_comparison must be a character string" = is.character(sample_comparison) & length(sample_comparison) == 1)

  # validate genes_to_display

  stopifnot("genes_to_display must be a character vector" = is.character(genes_to_display) | is.null(genes_to_display))


  # prepare volcano plots

  # filter by sample_comparison
  volcano_data <- x$dsr[[sample_comparison]]


  # # fix gene name

  # if(repair_genes){
  #   volcano_data <- volcano_data %>%
  #     purrr::map( function(x){
  #     gn <- rownames(x)
  #     gsym <- stringr::str_replace(gn, pattern = ".*:", replacement = "")
  #     gid <- stringr::str_replace(gn, pattern = ":.*", replacement = "")
  #     gn_new <- dplyr::if_else(duplicated(gsym) , gid ,
  #                              gsym)
  #     rownames(x) <- gn_new
  #     return(x)
  #   })
  # }
  #

  # generate plot

  pp <- parcutils::EnhancedVolcano2(toptable = volcano_data,
                                    x = "log2FoldChange" ,
                                    y = "pvalue",
                                    lab = rownames(volcano_data),
                                    selectLab = genes_to_display,
                                    FCcutoff = log2fc_cutoff,
                                    pCutoff = pval_cutoff,
                                    labSize = lab_size,
                                    pointSize =point_size,
                                    col_by_regul = T,
                                    col_up = col_up,
                                    col_down = col_down,
                                    col_other = col_other,
                                    title = sample_comparison , ...)

  return(pp)


}




#' Group replicates by samples.
#'
#' @param x an abject of class "parcutils". This is an output of the function [parcutils::run_deseq_analysis()].
#'
#' @return a tbl.
#' @export
#' @keywords internal.
#' @examples
#' \dontrun{
#' // TO DO
#' }
group_replicates_by_sample <- function(x){

  parcutils:::validata_parcutils_obj(x)

  sample_info <- purrr::map(parcutils:::get_all_named_expression_matrix(x), ~ ..1 %>% colnames() %>% .[-1]) %>%
    tibble::tibble(groups = names(.) , samples = .) %>%
    tidyr::unnest(cols = "samples")

  return(sample_info)
}





#' Visualise alignemnt stats in a bar plot
#'
#' @param x a named character vector denoting paths to *Log.final.out file generated by STAR aligner.
#' Names will be used to label Y axis in the plot. If not named, file names will be used to label.
#' @param col_total_reads a character string denoting bar color total reads. Default "#007e2f".
#' @param col_mapped_reads a character string denoting bar color total reads. Default "#ffcd12".
#' @param is_paired_data  a logical, default TRUE, denoting weather data is paired or single end.
#' By default STAR counts a paired-end read as one read. When `is_paired_data` set to `TRUE` it
#' doubles the count of total and mapped reads.
#'
#' @return a bar plot.
#' @export
#'
#' @examples
#'
#' file_dir <- system.file("extdata",package = "parcutils")
#'
#' file_paths <- fs::dir_ls(path = file_dir,
#' glob = "*Log.final.out" , recurse = T, type = "file") %>% as.character()
#'
#' get_star_align_log_summary_plot(x  = file_paths)
#' get_star_align_log_summary_plot(x  = file_paths, col_total_reads = "red",
#' col_mapped_reads = "blue")
#'
#'
get_star_align_log_summary_plot <- function(x,is_paired_data = TRUE,
                                            col_total_reads = "#007e2f" ,
                                            col_mapped_reads = "#ffcd12") {

  if(! (fs::file_exists(x) %>% all()) ){
    not_exist_indx <- which(fs::file_exists(x) == F)
    not_exist_file <- paste0(x[not_exist_indx], collapse = ",")
    stop(glue::glue("file(s) {not_exist_file } not exist."))
  }

  # assign names from path.
  if(is.null(names(x))) {
    names(x) <- stringr::str_replace(string = x, pattern = ".*/","")
  }

  star_summary <- x %>%
    purrr::map(~ parcutils::get_star_align_log_summary(..1)) %>%
    purrr::map_df(~ dplyr::bind_rows(..1), .id = "sample") %>%
    tidyr::pivot_wider(names_from = "type", values_from = "val")

  align_sumary <- star_summary  %>% dplyr::select(1,2,3,4)

  ##  format data for plot
  align_sumary_long <- align_sumary %>%

    # convert char to numeric where possible
    readr::type_convert() %>%

    # rename column names
    dplyr::rename(`Total reads` = `Number of input reads` ,
                  `Mapped reads` = `Uniquely mapped reads number`) %>%

    # order samples
    dplyr::mutate(sample = factor(sample, levels = sample %>% rev())) %>%

    # wide to long
    tidyr::pivot_longer(cols  = c(`Total reads`,`Mapped reads`),
                        names_to = "read_status" , values_to = "number_of_reads")

  if(is_paired_data){
    align_sumary_long <- align_sumary_long %>%
      dplyr::mutate(number_of_reads = number_of_reads * 2)
  }

  ## plot

  gp_align_summary <- align_sumary_long %>%
    ggplot2::ggplot(ggplot2::aes(x = sample , y =  number_of_reads , fill = read_status )) +
    ggplot2::geom_bar(stat = "identity" , position = "dodge") +
    ggplot2::scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6)) +
    ggplot2::scale_fill_manual(breaks = c("Total reads","Mapped reads"),
                               values = c(col_total_reads,col_mapped_reads) ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() + ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2:: xlab("Samples")

  if(is_paired_data){
    gp_align_summary <- gp_align_summary +
      ggplot2::ylab("Total Number of Reads (R1 + R2)")
  } else{
    gp_align_summary <- gp_align_summary +
      ggplot2::ylab("Total Number of Reads")
  }
  return(gp_align_summary)


}


# plot diff genes counts

#' Generate a barplot of DEG counts.
#'
#' @param x an object of class parcutils
#' @param col_up a character string, default #a40000, denoting valid a color name for "Up" regulated genes.
#' @param col_down a character string, default #16317d, denoting valid a color name "Down" regulated genes.
#' @param font_size a numeric, default 12, denoting font size in the plot.
#'
#' @return a ggplot
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
#'                          group_denominator = c("control"),
#'                          column_samples = c("control_rep1", "treat1_rep1", "treat2_rep1", "control_rep2", "treat1_rep2", "treat2_rep2", "control_rep3", "treat1_rep3", "treat2_rep3"))
#'
#'get_diff_gene_count_barplot(res)
get_diff_gene_count_barplot <- function(x,
                                        col_up="#a40000",
                                        col_down="#16317d",
                                        font_size = 12){

  validata_parcutils_obj(x)

  gp <- x$deg_summmary %>%
    tibble::enframe(name = "comparison" ,
                    value = "deg_count") %>%
    tidyr::unnest(cols = "deg_count") %>%
    dplyr::filter(regul != "other") %>%
    dplyr::mutate(regul = forcats::fct_relevel(regul, c("Up","Down")))
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

  return(gp)

}


#' Generate a correlation heat plot.
#' @description For all (or selected) samples this function generate a correlation heat box type plot.
#' @param x an object of class parcutils.
#' @param samples a character vector denoting samples to plot in scatter plot, default \code{NULL}. If set to NULL all samples are accounted.
#' @param genes a character vector denoting genes to consider in scatter plot, default \code{NULL}. If set to NULL all genes are accounted.
#' @param corr_method a character string, default \code{"pearson"}, denoting a value for correlation method. Value can be one of these \code{"pearson", "kendall", "spearman"}.
#' @param vis_method a character string, default \code{"full"}, denoting a value type of visualization. Value can be one of these \code{"square", "circle"}.
#' @param plot_type a character string, default \code{"square"}, denoting value for plot type. Value can be one of these \code{"full", "lower", "upper"}.
#' @param cluster_samples a logical, default \code{TRUE}, denoting weather to cluster samples or not.
#' @param show_diagonal a logical, default \code{TRUE}, denoting weather to show diagonal values or not.
#' @param show_corr_values a logical, default \code{TRUE}, denoting weather to show corr values or not.
#' @param col_corr_values  a character string, default \code{"yellow"}, denoting a valid color string for corr values.
#' @param size_corr_values a numeric, default \code{5}, denoting a size for corr values.
#' @param scale_range a numeric vector of length two denoting minimum and maximum value for the color scale, default \code{NULL}.
#'
#' @return a corr plot.
#' @export
#' @import ggeasy scales ggcorrplot viridis
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
#'                          group_denominator = c("control"),
#'                          column_samples = c("control_rep1", "treat1_rep1", "treat2_rep1", "control_rep2", "treat1_rep2", "treat2_rep2", "control_rep3", "treat1_rep3", "treat2_rep3"))
#'
#'get_corr_heatbox(res,samples = c("treatment1","control"),cluster_samples = FALSE,show_corr_values =T,
#'size_corr_values = 4)
get_corr_heatbox <- function(x,
                             samples = NULL,
                             genes = NULL,
                             corr_method = "pearson" ,# "pearson", "kendall", "spearman"
                             plot_type  = "full", # "full", "lower", "upper"
                             vis_method = "square", # square or circle
                             show_diagonal = TRUE,
                             show_corr_values = FALSE,
                             col_corr_values = "yellow",
                             size_corr_values = 5,
                             cluster_samples=FALSE,
                             scale_range = NULL

){

  # validate x
  parcutils:::validata_parcutils_obj(x)

  # validate samples
  stopifnot("'samples' must be a character vector or NULL." = is.character(samples) | is.null(samples))

  # validate genes
  stopifnot("'genes' must be a character vector or NULL." = is.character(genes) | is.null(genes))

  # validate corr_method
  match.arg(corr_method , c("pearson", "kendall", "spearman"))

  # validate vis method
  match.arg(vis_method , c("square","circle"))

  # validate plot type

  match.arg(plot_type , c("full","lower" ,"upper"))

  # validate cluster_samples

  stopifnot("'cluster_samples' must be a logical value. " = is.logical(cluster_samples))

  # validate show_diagonal
  stopifnot("'show_diagonal' must be a logical value." = is.logical(show_diagonal))

  # validate show_corr_values
  stopifnot("'show_corr_values' must be a logical value." = is.logical(show_corr_values) )

  # validate col_corr_values
  stopifnot("'col_corr_values must be a character string denoting a valid color name.'" =
              is.character(col_corr_values) & length(col_corr_values)==1)

  # validate size_corr_values
  stopifnot("'size_corr_values' must be a numeric value." = is.numeric(size_corr_values) &
              length(size_corr_values) == 1)

  #  validate scale range
  stopifnot("'scale_range' must be NULL or a numeric vector of length two."=
              is.null(scale_range) | (is.numeric(scale_range) & length(scale_range) == 2))

  # get expression values
  norm_expr_mat <- parcutils::get_normalised_expression_matrix(x = x,
                                                               samples = samples,
                                                               genes = genes,
                                                               summarise_replicates = F)

  # convert log2
  if(TRUE){
    norm_expr_mat <- norm_expr_mat %>% dplyr::mutate_if(is.numeric, ~ log2(. + 1))
  }

  # convert numeric df
  gene_id_column  <- colnames(norm_expr_mat)[1]
  norm_expr_mat <- norm_expr_mat %>% tibble::column_to_rownames(gene_id_column)

  corr_mat <- stats::cor(norm_expr_mat,method = corr_method)

  chb <- ggcorrplot::ggcorrplot(corr =  corr_mat,
                                method = vis_method ,
                                type =  plot_type,
                                hc.order = cluster_samples,
                                show.diag =  show_diagonal,
                                lab_col = col_corr_values ,
                                lab_size = size_corr_values,
                                lab = show_corr_values
  )


  # fix scale range if it is NULL
  if(is.null(scale_range)){
    scale_min = min(corr_mat)
    scale_max = max(corr_mat)
    scale_range = c(scale_min, scale_max)
  }

  # color plot
  suppressMessages(
    chb <- chb +
      ggplot2::scale_fill_gradientn(colours = viridis::viridis(n = 10) %>% rev(),
                                    limits = scale_range , oob=scales::squish

      ) +
      ggeasy::easy_add_legend_title("Corr")

  )




  return(chb)

}






