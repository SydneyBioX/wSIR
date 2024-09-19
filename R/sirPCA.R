#' sirPCA
#'
#' @description
#' This function performs eigendecomposition on `(X^H)^t %*% W %*% (X^H)`,
#' where `X^H` is matrix of scaled slice means.
#'
#' @param sliced_data matrix of scaled slice means
#' @param maxDirections integer (default 50), number of directions we want in
#' our final low-dimensional Z.
#' @param W matrix of slice weights. Output of createWeightMatrix function.
#' @param varThreshold numeric proportion of eigenvalues of variance in
#' `t(X_H) %*% W %*% X_H` to retain.
#' Default is 0.99. Select higher threshold to include more dimensions,
#' lower threshold to include less dimensions.
#'
#' @return list containing: 1) "eigenvectors" matrix of eigenvectors of
#' `(X^H)^t %*% W %*% (X^H)`
#' 2) "d" integer, selected number of dimensions
#' 3) "evalues" vector of eigenvalues of `(X^H)^t %*% W %*% (X^H)`
#'
#' @keywords internal

sirPCA <- function(sliced_data,
                   maxDirections = 50,
                   W = diag(nrow(sliced_data)),
                   varThreshold = 0.99) {

  nslices <- nrow(sliced_data)
  sliced_data <- as.matrix(sliced_data)
  m1 <- .matMultArma(t(sliced_data), W)
  m <- .matMultArma(m1, sliced_data)

  eig_m <- .fastEigen(m)
  all_pc <- eig_m$vectors[, ncol(m):1]
  eig_m_values <- base::pmax(rev(eig_m$values), 0)

  propvariance_explained <- base::cumsum(eig_m_values)/sum(eig_m_values)
  d <- which(propvariance_explained >= varThreshold)[1]
  d <- min(d, maxDirections)
  return(list(evectors = all_pc[,seq_len(d)],
              d = d,
              evalues = eig_m$values))
}
