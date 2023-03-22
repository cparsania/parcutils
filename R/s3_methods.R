#' S3 print for object `parcutils`
#'
#' @param object
#'
#' @return
#' @export
# print <- function(object) {
#   UseMethod("print")
# }
#' @export print.parcutils
#' @export
#' @keywords internal
print.parcutils <- function(object) {
  if(TRUE){
    cli::cat_boxx(cli::style_bold("Summary of DE analysis"))
    cli::cat_line()
    cli::cli_text("Total number of comparison{?s}: {length(object$de_comparisons)} \n")
    cli::cli_text("Total number of genes used for DE analysis: {nrow(object$norm_counts[[1]][[1]])} \n")

    cli::cat_line()
    cli::cli_text(cli::style_bold(cli::col_br_cyan("DE gene counts by comparison...")))
    cli::cat_line()
    for(i in object$de_comparisons){
      #i <- object$de_comparisons[[1]]
      n_up_genes <- parcutils::get_genes_by_regulation(x = object,sample_comparisons = i, regulation = "up") %>% length()
      n_down_genes <- parcutils::get_genes_by_regulation(x = object,sample_comparisons = i, regulation = "down") %>% length()
      cat(cli::style_italic(i) ,"\n")
      cli::cat_bullet(cli::col_red(sprintf("number of up genes   : %d.", n_up_genes)))
      cli::cat_bullet(cli::col_green(sprintf("number of down genes : %d.", n_down_genes)))
      cli::cat_rule(width = 30)
      cli::cat_line()
    }

  }

}




#' S3 print for object `parcutils_ir`
#'
#' @param object
#'
#' @return
#' @export
# print <- function(object) {
#   UseMethod("print")
# }
#' @export print.parcutils_ir
#' @export
#' @keywords internal
print.parcutils_ir <- function(object) {
  if(TRUE){
    cli::cat_line()
    cat(sprintf("Total number of RI used for DE analysis are %d. \n", nrow(object$norm_counts[[1]][[1]])))
    cat(sprintf("Total number of comparisons are %d. \n", length(object$de_comparisons)))
    cli::cat_line()
    cat("Number of DE RI in each comparison:\n")
    cli::cat_line()
    for(i in object$de_comparisons){
      #i <- object$de_comparisons[[1]]
      n_up_genes <- parcutils::get_genes_by_regulation(x = object,sample_comparisons = i, regulation = "up") %>% length()
      n_down_genes <- parcutils::get_genes_by_regulation(x = object,sample_comparisons = i, regulation = "down") %>% length()
      cat(cli::style_bold(i) ,"\n")
      cli::cat_bullet(cli::col_red(sprintf("number of up RI   : %d.", n_up_genes)))
      cli::cat_bullet(cli::col_green(sprintf("number of down RI : %d.", n_down_genes)))
      cli::cat_rule(width = 30)
    }
  }

}


#' S3 print for object `parcutils_ase`
#'
#' @param object
#'
#' @return
#' @export
# print <- function(object) {
#   UseMethod("print")
# }
#' @export print.parcutils_ase
#' @export
#' @keywords internal
print.parcutils_ase <- function(object) {
  if(TRUE){
    cli::cat_boxx(cli::style_bold("Summary of ASE diff. analysis"))
    cli::cat_line()
    cli::cli_text("Total number of comparison{?s}: {length(object$de_comparisons)} \n")
    cli::cli_text("Total number of ASE used for diff. analysis: {nrow(object$avg_psi[[1]][[1]])} \n")

    cli::cat_line()
    cli::cli_text(cli::style_bold(cli::col_br_cyan("Diff. ASE counts by comparison...")))
    cli::cat_line()
    for(i in object$de_comparisons){

      n_up_ase <- parcutils::get_ASE_by_regulation(x = object,sample_comparisons = i, regulation = "up") %>% length()
      n_down_ase <- parcutils::get_ASE_by_regulation(x = object,sample_comparisons = i, regulation = "down") %>% length()
      cat(cli::style_italic(i) ,"\n")
      cli::cat_bullet(cli::col_red(sprintf("number of up genes   : %d.", n_up_ase)))
      cli::cat_bullet(cli::col_green(sprintf("number of down genes : %d.", n_down_ase)))
      cli::cat_rule(width = 30)
      cli::cat_line()
    }

  }

}


