
#' Construct a parcutils_chip object
#' @description This function creates an object of the  {.cls parcutils_chip.} parcutils_chip is an extension of the the {.cls RangedSummarizedExperiment}.
#' This object make sure that it contains necessary elements to perform downstream chipseq analysis. Downstream exploratory and visualization functions will be based on this object.
#' @param x a list of the object(s) of {.cls normalisedMatrix}.
#' @param row_ranges an object of the class GenomicRanges. Each range must be identified by an unique id stored in a column \code{name},
#' @details RPM is calculated by the function RowSums applied on the each row of the
#' {.arg normalized_matrix}.
#' @return an object of the class parcutils_chip. Class parcutils_chip extends to RangedSummarizedExperiment. It contains data in the below format
#' \itemize{
#'  \item{"assay"} {These is the object of NormalisedMatrix for the sample.}
#'  \item{"columns"} {Columns are the columns in each assays. In this case they are the bins described in NormalisedMatrix.}
#'  \item{"rowData"} {RPM value for each feature (row) across samples. }
#' }
#' @export
#'
#' @examples
#'
#' signal = GRanges(seqnames = "chr1",
#' ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
#'                  end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)),
#' score = c(1, 2, 3, 1, 2, 3, 1, 2, 3))
#' target = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20), name = "feat_1")
#' x = normalizeToMatrix(signal, target, extend = 10, w = 2)
#'  rownames(x) <- "feat_1"
#' make_parcutils_chip(x, target)
#'
make_parcutils_chip <- function(x,row_ranges){

  # validate assays
  x <- .validate_a_list_of_normalised_matrix(x = x)

  # validate row_ranges
  meta_data_cols <- colnames(GenomicRanges::mcols(row_ranges))

  # check presence of the column `name`
  if(!any(meta_data_cols %in% "name")){
    cli::cli_abort("{.arg row_ranges} must contains feature id specified as a column {.val name}")
  }

  # check if feature ids are uniq

  if(any(duplicated(row_ranges$name))){
    cli::cli_abort(" feature id specified as a column {.val name} must have unique values.")
  }

  # rownames of each element (objects of normalisedMatrix) in assays and the `name` column of row_ranges must be identical

  rownames_equals_name  <- purrr::map_lgl(x, ~all( rownames(..1) == row_ranges$name))
  if(!all(rownames_equals_name)){
    rownames_differ <- names(rownames_equals_name)[!rownames_equals_name]
    cli::cli_abort("{.val rownames(x)} must be identical to the values in the column {.val name} of the {.arg row_ranges}")
  }


  # get RPM values

  feature_rpm <- normalised_matrix_to_rpm(x = x)

  # append RPM values as columns to original row_ranges object.

  row_ranges_rpm <- row_ranges %>%
    as.data.frame() %>%
    dplyr::left_join(feature_rpm, by = c("name" = "feature_id")) %>%
    plyranges::as_granges()

  # make se
  se <- SummarizedExperiment::SummarizedExperiment(assays = x,
                                                   rowData = row_ranges_rpm)
  se <- as(se,"parcutils_chip")
  return(se)
}


#' Set class parcutils_chip
#'
#' @return an object of class parcutils_chip
#' @export
#' @keywords internal
#' @examples
#' \dontrun{
#' // TO DO
#' }
parcutils_chip <- setClass(Class = "parcutils_chip", contains="RangedSummarizedExperiment")



#' Convert normalised matrix into RPM signal
#'
#' @param x a list of the object(s) of {.class normalisedMatrix}.
#' @details RPM is calculated by the function RowSums applied on the each row of the
#' {.arg x}.
#' @return a dataframe of RPM values. Rows are the features and columns are samples. Row names will be set from the rownames of the {.arg x}. Column names are the attribute {.code names} the list {.arg x}.
#' @export
#'
#' @examples
#' \dontrun{
#' // TODO
#' }
normalised_matrix_to_rpm <- function(x){

  norm_mats <- .validate_a_list_of_normalised_matrix(x = x)

  # prepare RPM for each peak

  df_rpm <- furrr::future_map(norm_mats , ~ ..1 %>% rowSums() %>% set_names(rownames(..1))) %>%
    # make it wide format
    tibble::enframe() %>%
    tidyr::unnest(cols = c(value)) %>%
    dplyr::mutate(feature_id = names(value)) %>%
    tidyr::pivot_wider(names_from = "name", values_from = "value", id_cols = c("feature_id")) %>%
    as.data.frame() %>%
    tibble::column_to_rownames("feature_id")


  # add prefix "RPM-" in the column names

  names(df_rpm) <- stringr::str_c("RPM-",names(df_rpm))

  # convert tibble

  if(TRUE){
    df_rpm %<>% tibble::rownames_to_column("feature_id") %>% tibble::as_tibble()
  }

  return(df_rpm)

}

# a function to validate a list of normalised matrix.

#' Validate a list of normalised matrix.
#'
#' @param x a list of the object(s) of {.class normalisedMatrix}.
#' @return a list of the object(s) of {.class normalisedMatrix}.
#' @export
#'
#' @examples
#' \dontrun{
#' // TODO
#' }
.validate_a_list_of_normalised_matrix <- function(x){

  norm_mats = x

  # validate input when it is a list
  if(is.list(norm_mats)){

    # check if a list is named

    if(!.is_a_list_named(lst = norm_mats)){
      cli::cli_abort(message = "{.arg x} must be an object of {.class normalizedMatrix} or a named list of {.class normalizedMatrix.}")
    }

    # check if all list elements are of the class normalizedMatrix

    if(!all(purrr::map_lgl(norm_mats, ~..1 %>% is("normalizedMatrix")))){
      cli::cli_abort(message = "{.arg x} must be an object of {.class normalizedMatrix} or a named list of {.class normalizedMatrix.}")
    }
  }

  # validate input when it is not a  list

  if(!is.list(norm_mats)){
    if(!is(norm_mats, "normalizedMatrix")){
      cli::cli_abort(message = "{.arg x} must be an object of {.class normalizedMatrix} or a named list of {.class normalizedMatrix.}")
    }

    # convert into a named list. This is the case when a single object of normalisedMatrix supplied. It will be named as "first"

    norm_mats <- list(norm_mats) %>% rlang::set_names("first")
  }

  # validate row names. Each element of the list must have row names assigned.
  contains_row_names <- purrr::map_lgl(norm_mats, ~ !is.null(rownames(..1)))
  if(!all(contains_row_names)){
    row_names_absent <- names(norm_mats)[which(!contains_row_names)]
    cli::cli_abort(message = "{.emph {row_names_absent}} sample{?s} {?does/do} not contain rownames. ")
  }

  return(norm_mats)
}

.is_a_list_named <- function(lst) length(lst) == sum(names(lst) != "",na.rm = TRUE)

# parallelly import bw files

#' Import several bw files using parallel processing.
#'
#' @param dir_path_to_bw_files a character string denoting a valid path to a parent directory in which several bw files are stored in a sample wise directory. Along with a file name (without extension)  name of the corresponding sample direcory will be used as name attribute.
#' @details This function uses [furrr::future_map()] to import the bw files. Following code is recommended to activate the parallel processing.
#'
#' \code {future::plan(future::multisession(), workers = 50)}
#' \code {.parallel_import_bw_files("path/do/bw/files)}
#' \code {future::plan(future::sequential())}
#'
#' Number of availabe threads/workers can be identified by  [future::availableWorkers()] and
#' [future::availableCores()].
#'
#' @return a list of GRanges object. Each GRanges object is an output of the function [rtracklayer::import.bw()].
#' @export
#'
#' @examples
#' \dontrun{
#' // TODO
#' }
.parallel_import_bw_files <- function(dir_path_to_bw_files){

  bw_signal_files <- fs::dir_ls(path = dir_path,regexp = "*.bw",recurse = T)

  # name of parent directory + file name (without extension) will be assigned as the `name` attribute.
  sep = "-"
  names(bw_signal_files) <- stringr::str_glue(
    "{fs::path_dir(bw_signal_files) %>% fs::path_file()}{sep}{fs::path_file(bw_signal_files) %>% fs::path_ext_remove()}" )

  if(names(bw_signal_files) %>% duplicated() %>% any()){
    cli::cli_abort("Attribute `name` assigned to each bw path must be unique.\nIn current setting name of the parent directory + file name (without extension) is assigned as the `name` attribute.")
  }


  # import bw files

  tictoc::tic()
  cli::cli_process_start(msg = "Improting bw files...")
  bw_signals <- furrr::future_map(bw_signal_files, ~ {
    rtracklayer::import.bw(..1)
  })
  cli::cli_process_done()
  tictoc::toc()

  return(bw_signals)
}



# read chipseq peak file AKA TF target regions.


#' Import top features (ranked by score) from a bed file.
#'
#' @param bed_feature_file a character string denoting a valid bed file.
#' @param topn a numeric value, default 5000, denoting number of top features to keep. Features are ranked by a column  {.code score}.
#' @param center logical, denoting whether to align each feature by central position. If TRUE, default, it will return a ranges of single nucleotide (width = 0) denting a central position of each feature.
#' @param ... Other arguments pass to the function [rtracklayer::import.bed()].
#' @return an object of the {.class GRanges} containing {.arg topn} features.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' // TO DO.
#' }
import_topn_bed_features <- function(bed_feature_file,topn = 5000, center = TRUE,...){

  topn_features_gr = rtracklayer::import(bed_feature_file,...)  %>%
    plyranges::arrange(-.$score) %>%
    tibble::as_tibble()  %>%
    dplyr::distinct() %>%
    dplyr::slice(1:{{ topn }}) %>%
    plyranges::as_granges()

  if(center){
    topn_features_gr <- plyranges::mutate(plyranges::anchor_center(topn_features_gr), width = 0)
  }

  return(topn_features_gr)

}



#' Create a HeatmapList from parcutils_chip.
#'
#' @param x an object of the class parcutils_chip.
#' @param cluster_targets_by_rpm logical, denoting weather to cluster targets by RPM value.
#' @param cluster_by NULL (default) or character vector denoting RPM columns to be used for clustering.
#' @param n_clust a numeric, default 2, denoting number of clusters for kmeans clustering.
#' @param cluster_rows logical, default FALSE, denoting weather to cluster rows using default clustering. This must be always false as kmeans is implimented.
#' @param heatmap_columns NULL (default) or character vector denoting columns to be displayed in the heatmap.
#' @param heatmap_column_title NULL (default) or character vector denoting title for each heatmap column.
#' @param heatmap_pos_line logical, default FALSE, internally passed to \code{pos_line} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param heatmap_pos_line_gp logical, default FALSE, internally passed to \code{pos_line_gp} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param heatmap_colors internally passed to \code{col} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' Default \code{circlize::colorRamp2(breaks = c(0,0.5,1),colors = c('#fee6ce','#e6550d','#E31A1C'))}.
#' @param heatmap_top_annotations internally passed to \code{top_annotations} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param heatmap_axis_name internally passed to \code{axis_name} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param heatmap_axis_name_rot internally passed to \code{axis_name_rot} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param heatmap_axis_name_gp internally passed to \code{axis_name_gp} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#' @param heatmap_border internally passed to \code{border} argument of [EnrichedHeatmap::EnrichedHeatmap()].
#'
#' @return a HeatmapList.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
#'
make_enriched_heatmap_list <- function(x,
                                       cluster_targets_by_rpm = TRUE,
                                       cluster_by = NULL,
                                       n_clust = 2,
                                       cluster_rows = FALSE,
                                       heatmap_columns = NULL,
                                       heatmap_column_title = NULL,

                                       heatmap_pos_line = FALSE,
                                       heatmap_pos_line_gp = grid::gpar(lty = 2),
                                       heatmap_colors = circlize::colorRamp2(breaks = c(0,0.5,1),
                                                                             colors = c('#fee6ce','#e6550d','#E31A1C')),
                                       heatmap_top_annotations = NULL,
                                       heatmap_axis_name= c("-3KB","Summit","+3KB"),
                                       heatmap_axis_name_rot = 0,
                                       heatmap_axis_name_gp = gpar(fontsize = 10),
                                       heatmap_border = TRUE){


  # validate x

  if(!is(x ,"parcutils_chip")){
    cli::cli_abort("x must be an object of {.cls parcutils_chip}")
  }


  # validate cluster_by

  rpm_sample_names <- rowData(x) %>% colnames()
  rpm_sample_names <- rpm_sample_names[grepl("^RPM-",rpm_sample_names)] # keep only those which starts with RPM

  if(cluster_targets_by_rpm){
    if(is.null(cluster_by)){
      cluster_by = rpm_sample_names # if NULL cluster using all samples.
    } else if(!all(cluster_by %in% rpm_sample_names)){
      cli::cli_abort("{.arg cluster_by} must be either {.cls NULL} or a {.cls character} vector of {.code colnames(rowData(x))} ")
    }
  }

  # validate heatmap_column_title

  if(!is.null(heatmap_column_title)){
    if(!length(heatmap_column_title) <= length(heatmap_columns)){
      cli::cli_abort("{.arg heatmap_column_title} must have length <= length of {.arg heatmap_columns} ")
    }
  }


  # validate heatmap_columns

  heatmap_column_names <- assays(x) %>% names()

  if(!is.null(heatmap_columns)){
    if(!all(heatmap_columns %in% heatmap_column_names)){
      cli::cli_abort("{.arg heatmap_columns} must be either {.cls NULL} or a {.cls character} vector of {.code assays(x) %>% names()} ")
    }
  }

  # perform clustering

  for_clust <- rowData(x)[cluster_by] %>% as.data.frame()

  # find peaks with NA rpm value
  na_peaks <- for_clust %>% rowSums() %>% is.na() %>%  which()

  # replace NA peak with value

  for_clust[na_peaks,] <- 0

  # prepare a list of enriched heatmap --> an output of draw function


  # kmeans clustering
  n_clust = n_clust
  set.seed(123)
  peak_clust <-paste0("cluster", kmeans(for_clust, centers = n_clust)$cluster)
  peak_clust <- tibble::tibble(name = rownames(for_clust), clust = peak_clust)



  # # heatmap top annotations
  # hm_top_anno <- ComplexHeatmap::HeatmapAnnotation(
  #   enriched = EnrichedHeatmap::anno_enriched(
  #   height = grid::unit(30,"mm"),
  #   gp = list(lwd= 3, col = "red", fontsize=15),value = "mean",
  #   axis_param = list(
  #     side = "left",
  #     facing = "inside"
  #   )))


  # plot heatmaps

  heatmaps <- purrr::map(heatmap_columns, ~EnrichedHeatmap::EnrichedHeatmap(mat = assays(x)[[..1]] ,
                                                                            name = ..1,
                                                                            col = heatmap_colors,
                                                                            column_title = ..1,
                                                                            axis_name = heatmap_axis_name,
                                                                            top_annotation = heatmap_top_annotations,
                                                                            cluster_rows = cluster_rows,
                                                                            pos_line = heatmap_pos_line,
                                                                            pos_line_gp = heatmap_pos_line_gp,
                                                                            axis_name_rot = heatmap_axis_name_rot,
                                                                            border = heatmap_border))



  # assign names
  names(heatmaps)  <- heatmap_columns


  # update heatmap column titles
  if(!is.null(heatmap_column_title)){

    for(i in seq_along(heatmap_column_title)){
      heatmaps[[i]]@matrix_color_mapping@name <- heatmap_column_title[i]
      heatmaps[[i]]@column_title <- heatmap_column_title[i]
    }

  }


  # make a list of ComplexHeatmap( an object of the class HeatmapList)
  hm_list = NULL

  for(i in heatmaps){
    hm_list = hm_list + i
  }


  hm_list_draw <- ComplexHeatmap::draw(hm_list, split = peak_clust$clust)

  return(hm_list_draw)

}


















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
#' @param ... Other arguments pass to  [EnrichedHeatmap::EnrichedHeatmap()].
#' @return an output of the function [EnrichedHeatmap::EnrichedHeatmap()].
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
#' @import ComplexHeatmap
#' @import EnrichedHeatmap
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



#' Get a seqlogo for flanking region across five prime end of the genomic range.
#' @description given an object of  `GenomicRanges` the function generates a seqlogo
#' of the region flanked around five prime end.
#' @param x an object of the class GenomicRanges
#' @param y an object of class BSgenome, default `BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38`
#' @param extend an integer, default 30, denoting number of basepairs to extend from five prime end of the range.
#'
#' @return a list containing ggseqlogo and sequences used to generate the logo.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
#'
get_five_prime_flank_motif <- function(x,
                                       y = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 ,
                                       extend = 30){

  if(!isClass(x,"GenomicRanges")){
    stop("x must a GenomicRanges object.")
  }

  stopifnot("'extend' must be an integer." = is.numeric(extend))


  # anchor 5 end and extend to the length 'l'
  extend_five_end <-  plyranges::mutate(plyranges::anchor_5p(x),width = extend)

  # anchor 3 end and extend to the length 'l'

  extend_three_end <- plyranges::stretch(plyranges::anchor_3p(extend_five_end),extend = extend)
  target_seq <- parcutils:::.map_granges_metadata(extend_three_end,
                                                  bs_genome_object = y) %>% .$seq

  #names(target_seq) <- extend_three_end$event_name

  for_motif <- list(as.character(target_seq))

  m1 <- ggplot2::ggplot() + ggseqlogo::geom_logo(for_motif) +
    ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks = 1:(extend*2), labels = c(-(extend):-1,1:extend))

  return(list(logo= m1, seq = for_motif))

}



#' Get a seqlogo for flanking region across three prime end of the genomic range.
#' @description given an object of  `GenomicRanges` the function generates a seqlogo
#' of the region flanked around three prime end.
#' @param x an object of the class GenomicRanges
#' @param y an object of class BSgenome, default `BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38`
#' @param extend an integer, default 30, denoting number of basepairs to extend from three prime end of the range.
#'
#' @return a list containing ggseqlogo and sequences used to generate the logo.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' }
#'
get_three_prime_flank_motif <-  function(x,
                                        y = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 ,
                                        extend = 30){

  if(!isClass(x,"GenomicRanges")){
    stop("x must a GenomicRanges object.")
  }

  stopifnot("'extend' must be an integer." = is.numeric(extend))

  # anchor 5 end and extend to the length 'l'
  extend_five_end <-  plyranges::mutate(plyranges::anchor_3p(x),width = extend)

  # anchor 3 end and extend to the length 'l'

  extend_three_end <- plyranges::stretch(plyranges::anchor_5p(extend_five_end),extend = extend)
  target_seq <- parcutils:::.map_granges_metadata(extend_three_end, bs_genome_object = y) %>% .$seq

  #names(target_seq) <- extend_three_end$event_name

  for_motif <- list(as.character(target_seq))

  m1 <- ggplot2::ggplot() + ggseqlogo::geom_logo(for_motif) +
    ggplot2::theme_bw() + ggplot2::scale_x_continuous(breaks = 1:(extend*2), labels = c(-(extend):-1,1:extend))

  return(list(logo= m1, seq = for_motif))

}
