
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


# vocano plot

#' Wrapper around EnhancedVolcano::EnhancedVolcano
#'
#' @param toptable \code{toptable} \link{EnhancedVolcano::EnhancedVolcano}
#' @param lab \code{lab} \link{EnhancedVolcano::EnhancedVolcano}
#' @param x \code{x} \link{EnhancedVolcano::EnhancedVolcano}
#' @param y \code{y} \link{EnhancedVolcano::EnhancedVolcano}
#' @param pCutoff \code{pCutoff} \link{EnhancedVolcano::EnhancedVolcano}
#' @param FCcutoff \code{FCcutoff} \link{EnhancedVolcano::EnhancedVolcano}
#' @param col_by_regul lgl, default T, whether to color variable by gene regulation - Up, Down, Other
#' @param ... Other parameters pass to \link{EnhancedVolcano::EnhancedVolcano}
#'
#' @return a volcano plot
#' @export
EnhancedVolcano2 <- function(toptable, lab, x, y, pCutoff = 10e-4,FCcutoff = 1.5,col_by_regul = T,...){

  if(col_by_regul) {

    keyvals <- ifelse(
      toptable[[x]] <= -FCcutoff  & toptable[[y]] <= pCutoff , '#2166ac',
      ifelse(toptable[[x]] >= FCcutoff &  toptable[[y]] <= pCutoff, "#b2182b",
             '#e0e0e0'))
    keyvals[is.na(keyvals)] <- '#e0e0e0'
    names(keyvals)[keyvals == '#2166ac'] <- 'Down'
    names(keyvals)[keyvals == '#e0e0e0'] <- 'Other'
    names(keyvals)[keyvals == '#b2182b'] <- 'Up'
  }

  plot  <- EnhancedVolcano::EnhancedVolcano(toptable = toptable, lab = lab,x = x,y = y,
                                            pCutoff = pCutoff, FCcutoff = FCcutoff, colCustom  = keyvals, ...)

  plot
}
