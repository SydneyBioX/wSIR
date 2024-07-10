#' subset_lower_tri
#'
#' @description
#' This function is within the code for finding correlation of distances
#'
#' @param m to fill
#'
#' @return to fill
#'
#' @keywords internal

subset_lower_tri = function(m) {
  mm = m[lower.tri(m, diag = FALSE)]
  return(c(mm))
}
