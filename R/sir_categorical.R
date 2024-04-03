# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param X no arguments
#' @param Y to fill
#' @param directions to fill
#' @param W to fill
#' @param ... to fill
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
sir_categorical <- function(X, Y, directions = 50, W = NULL, ...) {
  # X is data matrix
  # Y is a data.frame with a column "coordinate" that should be split
  # directions number of directions
  # W optionally provided and is output of cells_weight_matrix2()

  # do the transformation
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

  pc_dirs <- sir_PCA(sliced_data, directions = directions, W = W, ...)

  betas <- backsolve(R, pc_dirs$evectors)
  betas <- apply(betas, 2, function(x) x/sqrt(sum(x^2)))
  final_XB <- as.matrix(X) %*% betas


  return(list(scores = final_XB,
              directions = betas,
              estd = pc_dirs[[2]],
              W = W,
              evalues = pc_dirs$evalues)) # 1 is the transformed X, 2 is the rotation (for new data), 3 is the number of selected dimensions, 4 is the W
}
