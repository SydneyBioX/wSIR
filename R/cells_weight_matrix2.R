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
cells_weight_matrix2 <- function(coords, labels, alpha = 1) {

  avg_coords <- slicer_categorical(coords, labels)

  dist_mat <- as.matrix(dist(avg_coords, diag = TRUE, upper = TRUE))

  dist_norm <- 1 - dist_mat / max(dist_mat, na.rm = TRUE)

  weight_mat <- (dist_norm^alpha * 2) - 1

  return(weight_mat)
}
