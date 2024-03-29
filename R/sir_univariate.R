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
sir_univariate <- function(X, Y, directions = 10, categorical = FALSE, slices = 10, alpha = 1, W) {
  sliced_data <- slicer(X = X, Y = Y, slices = slices, categorical = categorical) # create sliced and averaged data

  ## sir_univariate

  ### inputs:
  # X: dataframe or matrix of dimensions n * p
  # Y: dataframe of size n * 1 (univariate response only), whose values are categorical or continuous (see "categorical" input)
  # directions: integer for number of SIR directions you wish to return
  # categorical: binary for if the response values are categorical or not
  # slices: number of slices to be used for SIR computation. Not needed if response is categorical, since in that case
  # number of slices = number of categories
  # alpha: integer, tuning parameter raising each entry of W to some power before the entries are scaled to [-1,1].
  # Motivation: further stretches away distant tiles, draws together similar ones.
  # Default value is 1, but best to adjust.
  # W: weight matrix to represent similarity between the slices. Only use for categorical response.
  # if none is provided, then (scaled) identity matrix is used. Argument alpha can still be used in that situation.

  ### outputs:
  # [[1]]: Z: low-dimensional representation of X. Dimensions n * d, d = directions argument.
  # [[2]]: B: SIR directions of size p * d (d = directions). To be used to project new data X_new of dimensions m * p into
  # low-dimensional space by performing X_new %*% B.


  # define weight matrix if not provided

  if (missing(W)) {
    if (categorical) {
      W = diag(table(Y[,1]))/nrow(Y)
      W = W^alpha
    } else {
      W = diag(slices)/nrow(X)
      W = W^alpha
    }
  }

  pc_dirs <- sir_PCA(sliced_data, directions = directions, W = W)

  betas <- multiplier(data = X, pc_dirs = pc_dirs[[1]])

  final_XB <- as.matrix(X) %*% betas

  return(list(final_XB, betas, pc_dirs[[2]])) # 1 is the transformed X, 2 is the rotation (for new data), 3 is the number of selected dimensions
}
