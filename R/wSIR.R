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
wSIR = function(X,
                coords,
                groups = rep(1, nrow(coords)),
                slices = 8,
                alpha = 4,
                maxDirections = 50,
                ...) {
  # ... passed onto sir_PCA e.g. argument varThreshold
  #browser()
  coords_split = split.data.frame(coords, groups)

  tile_allocations = lapply(coords_split, spatial_allocator2, slices = slices)

  tile_allocation = do.call(rbind, tile_allocations)

  tile_allocation$coordinate <- paste0(tile_allocation$coordinate,
                                       ", ",
                                       as.integer(factor(groups)))

  # tile_allocation <- spatial_allocator2(coords = coords, slices = slices)

  sliceName = "coordinate"
  labels = tile_allocation[,sliceName,drop = FALSE]
  
  
  # This Dmatrix is for the proportion
  H = table(tile_allocation$coordinate)
  Dmatrix <- diag(sqrt(H)/nrow(X), ncol = length(H))
  
  if (weighted) {
    corrMatrix = cells_weight_matrix2(coords, labels = labels, alpha = alpha)
  } else {
    corrMatrix = diag(ncol = length(H), nrow = length(H))
  }
  
  W <- Dmatrix %*% corrMatrix %*% Dmatrix
  
  
  wsir_obj <- sir_categorical(X = X,
                              Y = tile_allocation,
                              directions = maxDirections,
                              W = W, ...)
  return(wsir_obj)
}
