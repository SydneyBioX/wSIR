#' wSIR_specified_params
#'
#' @description
#' wSIR function for specified parameters alpha, slices. To be ued internally in wSIROptimisation and wSIR functions, each for
#' specific values of alpha and slices.
#'
#' @param X matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param group a factor indicating group level information for cells within and across samples (e.g. cell type).
#' @param samples sample ID of each cell. In total, must have length equal to the number of cells. For example, if
#' your dataset has 10000 cells, the first 5000 from sample 1 and the remaining 5000 from sample 2, you would write
#' samples = c(rep(1, 5000), rep(2, 5000)) to specify that the first 5000 cells are sample 1 and the remaining are sample 2.
#' Default is that all cells are from sample 1. Sample IDs can be of any format: for the previous example, you could write
#' samples = c(rep("sample 1", 5000), rep("sample 2", 5000)), and the result would be the same.
#' @param slices integer for number of slices on each axis of tissue. For example, slices = 4 creates 4 slices along
#' each spatial axis, yielding 4^2 = 16 tiles. Default is 8, suggested minimum of 3. Suggest to tune this parameter.
#' @param alpha integer to specify desired strength of spatial correlation. Suggest to tune this parameter on testing
#' dataset among e.g values c(0,2,4,8). alpha = 0 gives SIR implementation. Larger values give stronger spatial correlations.
#' @param ... arguments passed to sirCategorical
#'
#' @return wSIR object which is a list containing 5 (named) positions. 1) scores matrix containing low-dimensional
#' representation of X from wSIR algorithm. Dimension n * d, where d is chosen to include at least varThreshold proportion of variance.
#' 2) directions matrix containing WSIR directions, dimension p * d. Use this to project new data into low-dimensional space via X_new %*% directions.
#' 3) estd integer stating how many dimensions in the computed low-dimensional space are needed to account for varThreshold
#' proportion of variance. Same as number of dimensions in scores.
#' 4) W matrix weight matrix from createWeightMatrix function
#' 5) evalues vector containing p eigenvalues of t(X_H) %*% W %*% X_H. varThreshold parameter works on these evalues,
#' such that e.g the first j directions are included if the sum of the first j evalues equals 0.95% of the sum of all evalues.
#'
#' @importFrom methods is
#' @keywords internal
wSIRSpecifiedParams = function(X,
                               coords,
                               group = NULL,
                               samples = rep(1, nrow(coords)),
                               slices = 8,
                               alpha = 4,
                               # maxDirections = 50,
                               # varThreshold = 0.95
                               ...
                               ) {

  if (!is.null(group)) {

    if (methods::is(group, "factor")) {
      message("group is not a factor, setting as.factor")
      group <- factor(group)
    }

    coords2 = cbind(coords, group = group)
  } else {
    coords2 = coords
  }

  coords_split = split.data.frame(coords2, samples)

  tile_allocations = lapply(coords_split, spatialAllocator, slices = slices)

  tile_allocation = do.call(rbind, tile_allocations)

  tile_allocation$coordinate <- paste0(tile_allocation$coordinate,
                                       ", ",
                                       as.integer(factor(samples)))

  # tile_allocation <- spatial_allocator2(coords = coords, slices = slices)

  sliceName = "coordinate"
  labels = tile_allocation[,sliceName,drop = FALSE]


  # This Dmatrix is for the proportion
  H = table(tile_allocation$coordinate)
  Dmatrix <- diag(sqrt(H)/nrow(X), ncol = length(H))

  corrMatrix = createWeightMatrix(coords, labels = labels, alpha = alpha)

  W <- Dmatrix %*% corrMatrix %*% Dmatrix

  wsir_obj <- sirCategorical(X = X,
                             Y = tile_allocation,
                             W = W,
                             # maxDirections = maxDirections,
                             # varThreshold = varThreshold
                             ...
                             )

  return(wsir_obj)
}
