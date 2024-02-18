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
cells_weight_matrix <- function(coords, labels, alpha = 1) {

  ## cells_weight_matrix

  ### inputs:
  # coords: n * 2 matrix of spatial coordinates containing one row for each observation.
  # labels: vector of length n containing tile allocation of each observation, as produced by cells_weight_matrix.
  # alpha: integer, tuning parameter raising each entry of W to some power before the entries are scaled to [-1,1].
  # Motivation: further stretches away distant tiles, draws together similar ones.
  # Default value is 1, but best to adjust.

  ### output: Weight matrix of size (s^2) * (s^2) where s is the number of slices in each direction as
  # defined in the spatial_allocator function. Represents physical similarity between all pairs of tiles.


  slices <- length(unique(labels[,1]))
  empty_df <- matrix(rep(0, slices^2), nrow = slices) %>% as.data.frame()

  avg_locations <- slicer(X = coords, Y = as.data.frame(labels[,1]), categorical = TRUE)
  colnames(avg_locations) <- c("x", "y")

  for (i in 1:slices) {
    for (j in 1:slices) {
      x_dist <- avg_locations$x[i] - avg_locations$x[j]
      y_dist <- avg_locations$y[i] - avg_locations$y[j]
      dist_pair <- sqrt(x_dist^2 + y_dist^2)
      empty_df[i,j] = dist_pair
      empty_df[j,i] = dist_pair
    }
  }
  dist_df <- 1 - empty_df / max(empty_df)
  weight_mat <- dist_df %>% as.matrix()
  powered <- weight_mat^alpha
  return(powered*2-1)
}
