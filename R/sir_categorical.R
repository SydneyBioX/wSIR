# Hello world function

#' Day hello
#'
#' @description
#' This function performs WSIR based on provided gene expression data, tile allocations, and weight matrix.
#'
#' @param X matrix of normalised gene expression data for n genes across p cells. 
#' @param Y dataframe with 1 column named "coordinate" that is the tile allocation for each cell. There should be up to slices^2 unique tile IDs in this column.
#' @param directions Integer to specify maximum number of directions to retain in thee low-dimensional embedding of the data. Use if you need at most a certain number for a downstream task.
#' @param W Weight matrix created by cells_weight_matrix. Entry (i,j) represents the spatial correlation level between tiles i and j. The diagonal values should be all 1. If not provided, SIR implementation will be used. 
#' @param varThreshold numeric value specifying the desired proportion of variance to retain. If e.g 95% (default), then number of directions used will be such that their eigenvalues add up to 95% of the sum of all eigenvalues.
#'
#' @return list of outputs witth 5 named slots. They are the same as the output of the wSIR function: this is the final step in the wSIR function.
#'
#' @examples
#' #hello()
#'
#' @export
sir_categorical <- function(X, Y, directions = 50, W = NULL, varThreshold = 0.95) {

  # do the transformation (QR scaling method)
  n <- nrow(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  qr.Xc <- qr(Xc)
  R <- qr.R(qr.Xc)
  Z <- qr.Q(qr.Xc) * sqrt(n)

  sliced_data <- slicer_categorical(X = Z, Y = Y)

  # define weight matrix if not provided
  # want to find a better alternative to below
  if (is.null(W)) {
    W <- diag(table(Y$coordinate), ncol = nrow(sliced_data))/nrow(Y)
  }

  pc_dirs <- sir_PCA(sliced_data, directions = directions, W = W, varThreshold = varThreshold)

  betas <- backsolve(R, pc_dirs$evectors)
  betas <- apply(betas, 2, function(x) x/sqrt(sum(x^2)))
  final_XB <- as.matrix(X) %*% betas


  return(list(scores = final_XB,
              directions = betas,
              estd = pc_dirs[[2]],
              W = W,
              evalues = pc_dirs$evalues)) # 1 is the transformed X, 2 is the rotation (for new data), 3 is the number of selected dimensions, 4 is the W
}
