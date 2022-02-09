
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
#' \dontrun{
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
#' }
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

#' Wrapper around EnhancedVolcano::EnhancedVolcano
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


