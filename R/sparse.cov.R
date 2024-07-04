#' A function to calculate covariance of a sparse matrix. 
#'
#' @description
#' This function finds the covariance of a sparse matrix.
#'
#' @param x sparse matrix
#'
#' @return covariance matrix of x
#'
#' @keywords internal
sparse.cov <- function(x){
  n <- nrow(x)
  m <- ncol(x)
  ii <- unique(x@i)+1 # rows with a non-zero element

  Ex <- colMeans(x)
  nozero <- as.vector(x[ii,]) - rep(Ex,each=length(ii))        # colmeans

  covmat <- ( crossprod(matrix(nozero,ncol=m)) +
                crossprod(t(Ex))*(n-length(ii))
  )/(n-1)
  # sdvec <- sqrt(diag(covmat))
  # covmat/crossprod(t(sdvec))
  return(covmat)
}
