#' spatialAllocator
#'
#' @description
#' This function allocates each cell to a tile based on a specified number of tiles.
#'
#' @param coords dataframe contains the spatial position of each cell. Column names c("x", "y'). Must include row
#' names as integer from 1 to nrow(coords). This is automatically included in the wSIR function prior to this function.
#' @param slices integer the number of slices along each of the two spatial axes. The number of tiles will be
#' slices^2 since there are two spatial axes. E.g set slices = 4 to use 4^2 = 16 tiles. Default value is 3.
#'
#' @return output by itself is a matrix containing slice belonging for each axis in long format. When used in
#' lapply(coords_split, spatialAllocator, slices) as in wSIR the output is a dataframe with each cell's tile
#' allocation in the "coordinate" column.
#'
#' @keywords internal

spatialAllocator <- function(coords, slices = 3) {

  sliced <- lapply(coords, function(x) as.integer(cut(rank(x), slices)))
  allocation <- as.factor(do.call(paste, c(sliced, sep = ", ")))

  newcoords <- cbind(coords, coordinate = allocation)

  return(newcoords)
}
