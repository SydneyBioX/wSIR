#' wSIR
#'
#' @description
#' A function to perform supervised dimension reduction of a gene expression matrix using coordinates dataframe as response value. Incorporates weighting mechanism into SIR algorithm to generate low-dimensional representation of the data and allow for projection of new single-cell gene expression data into low-dimensional space. 
#'
#' @param X matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y"). 
#' @param groups unsure purpose, to fill
#' @param slices integer for number of slices on each axis of tissue. For example, slices = 4 creates 4 slices along each spatial axis, yielding 4^2 = 16 tiles. Default is 8, suggested minimum of 3. Suggest to tune this parameter.
#' @param alpha integer to specify desired strength of spatial correlation. Suggest to tune this parameter on testing dataset among e.g values c(0,2,4,8). alpha = 0 gives SIR implementation. Larger values give stronger spatial correlations.
#' @param maxDirections integer for upper limit on number of directions to include in low-dimensional space. Use if you need less than a certain number of dimensions for downstream analyes
#' @param varThreshold numeric proportion of eigenvalues of variance in t(X_H) %*% W %*% X_H to retain. Default is 0.95. Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#'
#' @return wSIR object which is a list containing 5 (named) positions. 1) scores matrix containing low-dimensional representation of X from wSIR algorithm. Dimension n * d, where d is chosen to include at least varThreshold proportion of variance.
#' 2) directions matrix containing WSIR directions, dimension p * d. Use this to project new data into low-dimensional space via X_new %*% directions. 
#' 3) estd integer stating how many dimensions in the computed low-dimensional space are needed to account for varThreshold proportion of variance. Same as number of dimensions in scores.
#' 4) W matrix weight matrix from cells_weight_matrix2 function
#' 5) evalues vector containing p eigenvalues of t(X_H) %*% W %*% X_H. varThreshold parameter works on these evalues, such that e.g the first j directions are included if the sum of the first j evalues equals 0.95% of the sum of all evalues.
#' 
#' @examples
#' # to fill with some simulated spatial data from a package
#'
#' @export
wSIR = function(X,
                coords,
                groups = rep(1, nrow(coords)),
                slices = 8,
                alpha = 4,
                maxDirections = 50,
                varThreshold = 0.95) {
  #browser()
  coords_split = split.data.frame(coords, groups)

  tile_allocations = lapply(coords_split, spatial_allocator2, slices = slices)

  tile_allocation = do.call(rbind, tile_allocations)

  tile_allocation$coordinate <- paste0(tile_allocation$coordinate,
                                       ", ",
                                       as.integer(factor(groups)))

  # tile_allocation <- spatial_allocator2(coords = coords, slices = slices)

  sliceName = "coordinate"
  labels = tile_allocation[,sliceName,drop = FALSE]
  
  
  # This Dmatrix is for the proportion
  H = table(tile_allocation$coordinate)
  Dmatrix <- diag(sqrt(H)/nrow(X), ncol = length(H))
  
  if (weighted) {
    corrMatrix = cells_weight_matrix2(coords, labels = labels, alpha = alpha)
  } else {
    corrMatrix = diag(ncol = length(H), nrow = length(H))
  }
  
  W <- Dmatrix %*% corrMatrix %*% Dmatrix
  
  
  wsir_obj <- sir_categorical(X = X,
                              Y = tile_allocation,
                              directions = maxDirections,
                              W = W,
                              varThreshold = varThreshold)
  return(wsir_obj)
}
