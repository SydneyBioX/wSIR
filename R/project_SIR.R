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
project_SIR = function(sir, newdata) {

  ## Projecting new data
  # sir is list object output from `sir_univariate`
  # newdata is a cells x features matrix of expression

  proj = newdata %*% sir[[2]]
  return(proj)
}
