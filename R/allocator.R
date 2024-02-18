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
allocator <- function(nrows, slices = 8) {

  ## allocator

  ### inputs:
  # nrows: number of rows in the X dataframe/matrix
  # slices: integer that is number of slices for our SIR algorithm

  ### output: vector of slice allocations of length n (e.g (1,1,2,2) )



  rem <- nrows%%slices
  sizes <- rep(0, slices)
  for (i in 1:slices) {
    if (i <= rem) {
      sizes[i] <- nrows%/%slices + 1
    } else {
      sizes[i] <- nrows%/%slices
    }
  }
  allocation <- c()
  for (i in 1:slices) {
    allocation <- allocation %>% append(rep(i, sizes[i]))
  }
  return(allocation)
}
