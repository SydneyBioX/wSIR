#' sirPCA
#'
#' @description
#' This function performs eigendecomposition on (X^H)^t %*% W %*% (X^H), where X^H is matrix of scaled slice means.
#'
#' @param sliced_data matrix of scaled slice means
#' @param directions integer, number of directions we want in our final low-dimensional Z.
#' @param W matrix of slice weights. Output of createWeightMatrix function.
#' @param varThreshold numeric proportion of eigenvalues of variance in t(X_H) %*% W %*% X_H to retain.
#' Default is 0.95. Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#'
#' @return list containing: 1) "eigenvectors" matrix of eigenvectors of (X^H)^t %*% W %*% (X^H)
#' 2) "d" integer, selected number of dimensions
#' 3) "evalues" vector of eigenvalues of (X^H)^t %*% W %*% (X^H)
#'
#' @keywords internal

sirPCA <- function(sliced_data,
                    directions,
                    W = diag(nrow(sliced_data)),
                    varThreshold = 0.95) {

  nslices <- nrow(sliced_data)
  sliced_data <- as.matrix(sliced_data)
  m <- (t(as.matrix(sliced_data)) %*% W %*% sliced_data)
  eig_m <- eigen(m)
  all_pc <- eig_m$vectors
  eig_m_values <- pmax(eig_m$values, 0)
  propvariance_explained <- cumsum(eig_m_values)/sum(eig_m_values)
  d <- which(propvariance_explained >= varThreshold)[1]
  d <- min(d, directions) # directions is a maximum number of directions
  return(list(evectors = all_pc[,1:d],
              d = d,
              evalues = eig_m$values))
}
