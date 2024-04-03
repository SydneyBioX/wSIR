# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param x no arguments
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
sparse.cov <- function(x){
  # function to calculate covariance of a sparse matrix
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
