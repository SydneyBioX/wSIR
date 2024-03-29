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
sir_categorical <- function(X, Y, directions = 50, W = NULL, ...) {
  # X is data matrix
  # Y is a data.frame with a column "coordinate" that should be split
  # directions number of directions
  # W optionally provided and is output of cells_weight_matrix2()

  sliced_data <- slicer_categorical(X = X, Y = Y)

  # define weight matrix if not provided
  # want to find a better alternative to below
  if (is.null(W)) {
    W <- diag(table(Y$coordinate))/nrow(Y)
  }

  pc_dirs <- sir_PCA(sliced_data, directions = directions, W = W, ...)

  betas <- multiplier(data = X, pc_dirs = pc_dirs[[1]])

  final_XB <- as.matrix(X) %*% betas

  return(list(final_XB, betas, pc_dirs[[2]], W)) # 1 is the transformed X, 2 is the rotation (for new data), 3 is the number of selected dimensions, 4 is the W
}
