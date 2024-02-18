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
slicer_categorical <- function(X, Y) {

  # this function specifically slices with categorical Y
  # Y must have a column "coordinate" accessible via $coordinate
  # created by shila

  avg_X = do.call(rbind,lapply(split.data.frame(X, Y$coordinate), colMeans))

  return(avg_X)

}
