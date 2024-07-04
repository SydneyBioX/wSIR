#' cells_weight_matrix2
#'
#' @description
#' A function to create the weight matrix given the location of the cells, tile allocations and desired spatial weighting strength. Weight matrix entries represent level of spatial correlation between all pairs of tiles.
#'
#' @param coords dataframe of dimension n * 2. Column names c("x", "y"). Spatial position of each cell.
#' @param labels dataframe of dimension n * 1, column name c("coordinate"). Tile allocation of each cell. This is automatically created in the wSIR function.
#' @param alpha numeric value giving strength of spatial weighting matrix. alpha = 0 returns identity matrix and equals SIR. Large alpha values tend all entries towards 1. Default is 4.
#'
#' @return matrix containing the weight value for all pairs of tiles. Each value is between 0 and 1, with 1 always on the diagonal.
#'
#' @keywords internal

cells_weight_matrix2 <- function(coords, labels, alpha = 4) {
  alpha = 4/alpha # alpha_old = 0 (function argument = 0) gives alpha_new = Inf (rather than undefined) which equals SIR

  #browser()
  avg_coords <- slicer_categorical(coords, labels)

  avg_coords_groups <- gsub(".*, ", "", rownames(avg_coords))

  dist_mat <- as.matrix(stats::dist(avg_coords, diag = TRUE, upper = TRUE))

  # mask out the slices that are far away from each other
  D2 = as.matrix(stats::dist(avg_coords_groups, method = "manhattan"))
  D2[D2 > 0] <- Inf

  dist_mat_new = dist_mat + D2

  dist_norm <- (1 - dist_mat_new / max(dist_mat, na.rm = TRUE))^alpha
  weight_mat <- dist_norm
  weight_mat[!is.finite(weight_mat)] <- 0 # turn -inf into 0
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
