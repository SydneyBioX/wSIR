#' subsetLowerTri
#'
#' @description
#' This function is within the code for finding correlation of distances. It takes a distance matrix as an input,
#' and returns the lower triangle. This is useful because a distance matrix is symmetric, so we can take the lower
#' triangle only to reduce computation time
#'
#' @param m distance matrix
#'
#' @return vector of elements in the lower triangle of the distance matrix (input m)
#'
#' @keywords internal

subsetLowerTri = function(m) {
  mm = m[lower.tri(m, diag = FALSE)]
  return(c(mm))
}
