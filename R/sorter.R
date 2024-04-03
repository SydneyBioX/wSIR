# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param X no arguments
#' @param Y to fill
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
sorter <- function(X, Y) {

  # Hidden functions

  ## sorter

  ### inputs:
  # X: dataframe or matrix
  # Y: dataframe or matrix with one column

  ### output:
  # dataframe with rows sorted by their Y values



  # turn X and Y into data.frames
  X <- X %>% as.data.frame()
  Y <- Y %>% as.data.frame()
  # combine X and Y
  data <- cbind(X, Y)
  # give newly created Y column a name
  colnames(data)[ncol(data)] <- "y_vals"
  # sort by Y
  data <- data %>% arrange(y_vals)
  # remove column of Y values
  data <- data[,-ncol(data)]
  # return sorted dataframe
  return(data)
}
