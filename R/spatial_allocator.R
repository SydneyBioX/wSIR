# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param coords no arguments
#' @param slices to fill
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
spatial_allocator <- function(coords, slices = 3) {

  ## spatial_allocator

  ### inputs:
  # coords: n * 2 dataframe of spatial locations, one row per observation (cell)
  # note: the column names of coords MUST be "x" and "y" (lowercase, no quotation marks).
  # slices: integer for number of slices in each direction. E.g if you use slices = 3 then you will get
  # 3 slices in each of the x and y directions, leading to 3 x 3 = 9 tiles in total.

  ### output: vector of length n specifying the tile allocation of each observation.


  coords$id <- c(1:nrow(coords))
  coords <- coords %>% arrange(x)
  nX <- nrow(coords)
  x_allocation <- allocator(nrows = nX, slices = slices)
  coords$x_slice <- x_allocation

  coords <- coords %>% arrange(y)
  y_allocation <- allocator(nrows = nX, slices = slices)
  coords$y_slice <- y_allocation

  coords$coordinate <- paste0(coords$x_slice, ", ", coords$y_slice)

  coords <- arrange(coords, id)

  coords <- coords %>% dplyr::select(-c("x_slice", "y_slice", "id"))
  return(coords)
}
