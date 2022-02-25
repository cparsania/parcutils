#' S3 print for object `parcutils`
#'
#' @param object
#'
#' @return
#' @export
#' @keywords internal
print <- function(object) {
  UseMethod("print")
}

#' @export
#' @keywords internal
print.parcutils <- function(object) {
  if(TRUE){
    cli::cat_line()
    cat(sprintf("Total number of genes used for DEG analysis are %d. \n", nrow(object$norm_counts[[1]][[1]])))
    cat(sprintf("Total number of comparisons are %d. \n", length(object$comp)))
    cli::cat_line()
    cat("Number of DEGs in each comparison:\n")
    cli::cat_line()
    for(i in object$comp){
      #i <- deg_all_comparisons$comp[[1]]
      n_up_genes <- parcutils::get_genes_by_regulation(x = object,sample_comparisons = i, regulation = "up") %>% length()
      n_down_genes <- parcutils::get_genes_by_regulation(x = object,sample_comparisons = i, regulation = "down") %>% length()
      cat(cli::style_bold(i) ,"\n")
      cli::cat_bullet(cli::col_red(sprintf("number of up genes   : %d.", n_up_genes)))
      cli::cat_bullet(cli::col_green(sprintf("number of down genes : %d.", n_down_genes)))
      cli::cat_rule(width = 30)
    }

  }

}
