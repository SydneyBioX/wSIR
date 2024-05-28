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
  #browser()
  avg_coords <- slicer_categorical(coords, labels)

  avg_coords_groups <- gsub(".*, ", "", rownames(avg_coords))

  dist_mat <- as.matrix(dist(avg_coords, diag = TRUE, upper = TRUE))

  # mask out the slices that are far away from each other
  D2 = as.matrix(dist(avg_coords_groups, method = "manhattan"))
  D2[D2 > 0] <- Inf

  dist_mat_new = dist_mat + D2

  dist_norm <- (1 - dist_mat_new / max(dist_mat, na.rm = TRUE))^alpha
  weight_mat <- dist_norm
  #weight_mat <- (dist_norm^alpha* 2) - 1
  # ensure it is psd
  eig <- eigen(weight_mat)
  k <- eig$values > 1e-8
  weight_mat <- eig$vectors[, k, drop = F] %*%
    diag(eig$values[k]) %*%
    t(eig$vectors[, k, drop = F])

  weight_mat[!is.finite(weight_mat)] <- 0 #(-1)

  return(weight_mat)
}
