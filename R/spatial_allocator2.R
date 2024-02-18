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
spatial_allocator2 <- function(coords, slices = 3) {

  # spatial_allocator2 <- function(coords, slices = 3) {
  #
  #   sliced = apply(coords, 2, function(x){
  #     as.integer(cut(rank(x), slices))
  #   }, simplify = FALSE)
  #
  #   allocation = droplevels(do.call(interaction, c(sliced, sep = ", ")))
  #
  #   newcoords <- cbind(coords, coordinate = allocation)
  #
  #   return(newcoords)
  # }

  sliced <- lapply(coords, function(x) as.integer(cut(rank(x), slices)))
  allocation <- as.factor(do.call(paste, c(sliced, sep = ", ")))

  newcoords <- cbind(coords, coordinate = allocation)

  return(newcoords)
}
