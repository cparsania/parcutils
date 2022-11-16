#' Generate a barplot showing counts of splice events
#'
#' @param x an object of class \code{NxtSE}
#'
#' @return ggplot.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
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





