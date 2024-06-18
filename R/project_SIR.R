# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param sir no arguments
#' @param newdata to fill
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
project_wSIR = function(wsir, newdata) {

  ## Projecting new data
  # wsir is list object output from `wSIR`
  # newdata is a cells x features matrix of expression

  proj = newdata %*% wsir[[2]]
  return(proj)
}
