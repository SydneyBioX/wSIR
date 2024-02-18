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
means <- function(dataset_one_slice) {

  ## means

  ### input: dataframe of all rows and columns within one slice

  ### output: vector of means (length = number of columns in input) for mean across all rows for each column in input.


  nc <- ncol(dataset_one_slice)
  vals <- c()
  means_array <- colMeans(dataset_one_slice)
  for (i in 1:nc) {
    vals <- vals %>% append(means_array[[i]])
  }
  return(vals)
}
