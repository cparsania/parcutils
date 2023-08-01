#' Generate a ChIP signal heatmap for a given sample.
#' #' @description This function plots the ChIP signal heatmap for the given set of target region. ChIP signals shown in the plot will be plotted as an output of [EnrichedHeatmap::normalizeToMatrix()].
#' @details To obtain the order of targets from the output heatmap the suggested approach is `hm = ComplexHeatmap::draw(hm); ComplexHeatmap::row_order(hm)`.
#' @param path_to_bw_file  a character string denoting a path to a bw file. String can be an absolute path or the valid URL pointing to bw file.
#' @param target_region an object of the class [GenomicRanges]. Singal will be plotted against this regions.
#' @param extend pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param background pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param mean_mode pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param w pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param smooth pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param target_ratio pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param name pass to the same argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param col pass to the same argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param column_title pass to the same argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param axis_name pass to the same argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param top_annotation pass to the same argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param ...
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' ##
#' }
#'
get_chip_signal_heatmap <- function(path_to_bw_file ,
                                    target_region,
                                    extend = 5000,
                                    background = NA,
                                    mean_mode ="w0",
                                    w = 50,
                                    smooth = TRUE,
                                    target_ratio =0,
                                    name =NULL,
                                    col = NULL,
                                    column_title = NULL,
                                    axis_name = NULL,
                                    top_annotation = NULL,...
){


  x <- path_to_bw_file
  y <- target_region

  # load bw file
  bw_signal  <- rtracklayer::import.bw(x)

  # match sequence levels between signal and target

  GenomeInfoDb::seqlevelsStyle(bw_signal) <- "UCSC"
  GenomeInfoDb::seqlevelsStyle(y) <- "UCSC"

  # prepare a signal matrix for target regions

  bw_mat <- EnrichedHeatmap::normalizeToMatrix(signal = bw_signal,
                                               target = y,
                                               value_column = "score",
                                               background = background,
                                               extend = extend,
                                               mean_mode = mean_mode,
                                               keep = c(0, 0.99),
                                               w = w,
                                               target_ratio =target_ratio,
                                               smooth = smooth)

  # plot heatmap
  hm <-  EnrichedHeatmap::EnrichedHeatmap(mat = bw_mat ,
                                          name = name,
                                          col = col,
                                          column_title = column_title,
                                          axis_name = axis_name,
                                          top_annotation = top_annotation,...)

  return(hm)
}




#' Generate a ChIP signal heatmap of treatment over control.
#' @description This function plots the ChIP signal heatmap for the given set of target region.
#' ChIP signals shown in the plot is the ratio of the treatment over control. Ratio will be calculated for each bin of the target region. Prior to calculate the ratio integer 1 will be added to both numerator and denominator to avoid the Inf and NaN from the final matrix.
#' @details To obtain the order of targets from the output heatmap the suggested approach is `hm = ComplexHeatmap::draw(hm); ComplexHeatmap::row_order(hm)`.
#' @param path_control_bw a character string denoting a path to a bw file. String can be an absolute path or the valid URL pointing to bw file. Singals from this file will be used as numerator to calculate the fold change.
#' @param path_treatment_bw a character string denoting a path to a bw file. String can be an absolute path or the valid URL pointing to bw file. Singals from this file will be used as denominator to calculate the fold change.
#' @param target_region an object of the class [GenomicRanges]. Singal will be plotted against this regions.
#' @param name_control a chanracter string denoting a name of the control sample. If NULL, the name will be derieved from the argument 'path_control_bw'.
#' @param name_treatment a chanracter string denoting a name of the treatment sample. If NULL, the name will be derieved from the argument 'path_treatment_bw'.
#' @param background pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param extend pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param mean_mode pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param keep pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param w pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param target_ratio pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param smooth pass to the same argument of [EnrichedHeatmap::normalizeToMatrix()].
#' @param ... other arguments pass to the function [EnrichedHeatmap::EnrichedHeatmap()].
#' @return an object of the class Heatmap from the package ComplexHeatmap.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
#'
get_chip_signal_over_control_heatmap <- function(path_control_bw,
                                                 path_treatment_bw,
                                                 target_region = y,
                                                 name_control = NULL,
                                                 name_treatment = NULL,
                                                 background = background,
                                                 extend = extend,
                                                 mean_mode = mean_mode,
                                                 keep = c(0, 0.99),
                                                 w = w,
                                                 target_ratio =target_ratio,
                                                 smooth = smooth,...){


  # set names for control and treatment if it is NULL

  if(is.null(name_control)){
    name_control = fs::path_file(path_control_bw) %>% fs::path_ext_remove()
  }

  if(is.null(name_treatment)){
    name_treatment = fs::path_file(path_treatment_bw) %>% fs::path_ext_remove()
  }

  # make a named list
  bw_files <- list(control = path_control_bw, treatment = path_treatment_bw)
  names(bw_files) <- c(name_control, name_treatment)

  # read bw files
  bw_signals  <- purrr::map(bw_files, ~ rtracklayer::import(..1))

  # match sequence levels between signal and target

  for(i in seq_along(bw_signals)) {
    GenomeInfoDb::seqlevelsStyle(bw_signals[[i]]) <- "UCSC"
  }

  GenomeInfoDb::seqlevelsStyle(MBD3_Scr_IFNb_summit) <- "UCSC"


  # prepare a signal matrix for target regions
  bw_mats <- purrr::map(bw_signals , ~ EnrichedHeatmap::normalizeToMatrix(signal = ..1,
                                                                          target = target_region,
                                                                          value_column = "score",
                                                                          background = NA,
                                                                          extend = 5000,
                                                                          mean_mode = "w0",
                                                                          keep = c(0, 0.99),
                                                                          w = 50,
                                                                          target_ratio =0,
                                                                          smooth = TRUE))


  # calculate fold change treatment / control

  bw_mat_FC <- ((bw_mats[[name_control]]+1) / (bw_mats[[name_treatment]]+1))


  # plot heatmap
  hm <- EnrichedHeatmap::EnrichedHeatmap(mat = bw_mat_FC,
                                         col = hm_col,
                                         column_title = stringr::str_c(name_treatment,"/",name_control),
                                         axis_name = axis_name,
                                         top_annotation = top_annotation,...)


  return(hm)

}
