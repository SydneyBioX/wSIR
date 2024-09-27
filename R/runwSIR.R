#' runwSIR
#'
#' @description
#' Perform wSIR on cells, based on the expression data and a reducedDim in a
#' SingleCellExperiment or SpatialExperiment object
#'
#' @param x A numeric matrix of normalised gene expression data where rows are
#' features and columns are cells. Alternatively, a SingleCellExperiment or
#' SpatialExperiment containing such a matrix
#' @param name string to specify the name to store the result in the reducedDims
#' of the output. Default is "wSIR"
#' @param scores_only logical whether only the wSIR scores should be calculated.
#' If FALSE additional information about the wSIR model will be stored in the
#' attributes of the object. Default FALSE.
#' @param ... arguments passing to `calculatewSIR`
#'
#' @return If `x` is matrix-like, a list containing wSIR scores, loadings, etc.
#' If `x` is a SingleCellExperiment or SpatialExperiment, the same object is
#' returned with an additional slot in `reducedDims(..., name)` corresponding
#' to the wSIR scores matrix. If `scores_only = FALSE`, then the attributes of
#' the wSIR scores contain the following elements:
#' - directions
#' - estd
#' - W
#' - evalues
#'
#' @importFrom methods is
#' @importFrom SingleCellExperiment reducedDim
#'
#' @examples
#' data(MouseData)
#' library(SingleCellExperiment)
#' library(SpatialExperiment)
#'
#' sce = SingleCellExperiment(assays = list(logcounts = t(sample1_exprs)),
#' reducedDims = list(spatial = sample1_coords))
#'
#' sce = runwSIR(x = sce, dimred = "spatial")
#'
#' spe = SpatialExperiment(assays = list(logcounts = t(sample1_exprs)),
#' spatialCoords = as.matrix(sample1_coords))
#'
#' spe = runwSIR(x = spe, spatialCoords = TRUE)
#'
#' @export
runwSIR <- function(x,
                    name = "wSIR",
                    scores_only = FALSE,
                    ...) {

  isMatLike <- methods::is(x, "matrix")

  wsir_obj <- calculatewSIR(x = x, ...)

  if (isMatLike) {

    return(wsir_obj)

  }

  dr <- wsir_obj$scores

  if (scores_only) {

    base::attr(dr, "directions") <- wsir_obj$directions
    base::attr(dr, "estd") <- wsir_obj$estd
    base::attr(dr, "W") <- wsir_obj$W
    base::attr(dr, "evalues") <- wsir_obj$evalues

  }

  SingleCellExperiment::reducedDim(x, name) <- dr

  return(x)
}
