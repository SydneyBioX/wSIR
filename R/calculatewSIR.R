#' calculatewSIR
#'
#' @description
#' Perform wSIR on cells, based on the expression data and a reducedDim in a
#' SingleCellExperiment or SpatialExperiment object
#'
#' @param x A numeric matrix of normalised gene expression data where rows are
#' features and columns are cells. Alternatively, a SingleCellExperiment or
#' SpatialExperiment containing such a matrix
#' @param assay.type if `x` is a SingleCellExperiment of SpatialExperiment then
#' this is the assay for which wSIR will be calculated. Default "logcounts".
#' @param dimred String or integer scalar specifying the dimensionality reduction
#' slot for which to use for the slicing mechanism. Ignored if `coords` given.
#' @param colData_columns character vector specifying the subset of colData
#' columns to be used for the wSIR slicing mechanism. Ignored if `coords` or
#' `dimred` given
#' @param spatialCoords logical indicating if spatialCoords should be used for
#' the wSIR slicing mechanism. Ignored if `coords`, `dimred`, or
#' `colData_columns` given, or if `x` is not a SpatialExperiment object.
#' @param ... arguments passing to `wSIR`
#'
#' @return A wSIR object
#'
#' @importFrom SummarizedExperiment assay assays colData
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom methods is as
#' @importFrom BiocGenerics t
#'
#' @examples
#' data(MouseData)
#' library(SingleCellExperiment)
#' sce = SingleCellExperiment(assays = list(logcounts = t(sample1_exprs)),
#' reducedDims = list(spatial = sample1_coords))
#'
#' obj = calculatewSIR(x = sce,
#'   dimred = "spatial")
#'
#' @export
calculatewSIR <- function(x,
                          assay.type = "logcounts",
                          dimred = NULL,
                          colData_columns = NULL,
                          spatialCoords = FALSE,
                          ...) {

  isMatLike <- methods::is(x, "matrix")

  if (isMatLike) {

    wsir_obj <- wSIR(X = x, ...)

    return(wsir_obj)

  }

  if (!assay.type %in% names(SummarizedExperiment::assays(x))) {
    stop("assay.type not within assays of x")
  }

  X <- SummarizedExperiment::assay(x, assay.type)

  if (!is.null(dimred)) {

    if (!dimred %in% names(SingleCellExperiment::reducedDims(x))) {
      stop("dimred not within reducedDims of x")
    }

    coords <- SingleCellExperiment::reducedDim(x, dimred)

  } else {

    if (!is.null(colData_columns)) {

      if (!all(colData_columns %in% colnames(SummarizedExperiment::colData(x)))) {
        stop("not all colData_columns names are in colnames(colData(x))")
      }

      coords <- SummarizedExperiment::colData(x)[,colData_columns, drop = FALSE]
      coords <- methods::as(coords, "data.frame")

    } else {

      if (!is.null(spatialCoords)) {

        coords <- SpatialExperiment::spatialCoords(x)
        coords <- base::as.data.frame(coords)
        colnames(coords)[seq_len(2)] <- c("x", "y")

      }

    }

  }

  wsir_obj <- wSIR(X = BiocGenerics::t(X), coords = coords, ...)

  return(wsir_obj)
}
