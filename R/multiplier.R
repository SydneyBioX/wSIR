# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param data no arguments
#' @param pc_dirs to fill
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
multiplier <- function(data, pc_dirs) {

  ## multiplier

  ### inputs:
  # data: original X data
  # pc_dirs: eigenvectors of (X^H)^t %*% W %*% (X^H) as produced by sir_PCA

  ### output: matrix SIR directions from X, where each column of output is an SIR direction.


  if ("matrix" %in% class(data)) {
    cov_mat <- cov(data)
  } else {
    cov_mat <- sparse.cov(data)
  }
  multiplied <- solve(cov_mat, pc_dirs) # note: this will give an error if p >= n
  return(multiplied)
}
