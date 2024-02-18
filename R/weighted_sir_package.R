# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param none no arguments
#'
#' @return prints hello world
#'
#' @examples
#' hello()
#'
#' @export
weighted_sir_package <- function(X, coords, slices = 8, directions = 10, W = diag(slices)/nrow(X), alpha = 1) {


  ## weighted_sir_package

  ### inputs:
  # X: dataframe or matrix of dimensions n * p
  # coords: n * 2 matrix of spatial coordinates containing one row for each observation.
  # slices: number of slices to be used for SIR computation. Not needed if response is categorical, since in that case
  # number of slices = number of categories
  # directions: integer for number of SIR directions you wish to return
  # W: weight matrix to represent similarity between the slices, of size size (s^2) * (s^2)
  # where s is the number of slices in each spatial direction (x, y). If not provided then it will be created
  # automatically by the cells_weight_matrix function.
  # alpha: integer, tuning parameter raising each entry of W to some power before the entries are scaled to [-1,1].
  # Motivation: further stretches away distant tiles, draws together similar ones.
  # Default value is 1.

  ### outputs:
  # [[1]]: Z: low-dimensional representation of X. Dimensions n * d, d = directions argument.
  # [[2]]: B: SIR directions of size p * d (d = directions). To be used to project new data X_new of dimensions m * p into
  # low-dimensional space by performing X_new %*% B.

  tile_allocation <- spatial_allocator(coords = coords, slices = slices)
  if (missing(W)) {
    W <- cells_weight_matrix(coords = coords, labels = tile_allocation, alpha = alpha)
  }
  wsir_obj <- sir_univariate(X = X,
                             Y = tile_allocation,
                             directions = directions,
                             categorical = TRUE,
                             slices = slices,
                             alpha = alpha,
                             W = W)
  return(wsir_obj)
}
