#' MouseGastrulationData
#'
#' Data set consists of spatial transcriptomics data from a mouse embryo.
#' Two samples: first has gene expression data (11026 cells by 351 genes) and spatial
#' coords on the two-dimensional spatial plane, second contains gene expression data
#' only (8425 cells by 351 genes).
#'
#'
#' @name MouseData
#' @aliases sample1_exprs sample1_coords sample2_exprs
#' @docType data
#' @format \code{sample1_exprs} has a row for each cell in sample 1 and a column for the
#' expression level of each gene. \code{sample1_coords} has a row for each cell in sample 1
#' and a column for its position in each of the two spatial axes. \code{sample2_exprs} has
#' a row for each cell in sample 2 and a column for the expression level of each gene.
#' \code{sample2_coords} has a row for each cell in sample 2 and a column for its position
#' in each of the two spatial axes. \code{sample3_exprs} has a row for each cell in sample 3
#' and a column for the expression level of each gene. \code{sample3_coords} has a row for
#' each cell in sample 3 and a column for its position in each of the two spatial axes.
#' @source Integration of spatial and single-cell transcriptomic data elucidates mouse
#' organogenesis, \emph{Nature Biotechnology}, 2022.  Webpage:
#' \url{https://www.nature.com/articles/s41587-021-01006-2}
#' @keywords datasets
NULL
