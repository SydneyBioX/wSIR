#' project_wSIR
#'
#' @description
#' function to project new gene expression data into low-dimensional space
#'
#' @param wsir wsir object that is usually the output of wSIR function. If you want to project new data into low-dim space following a 
#' different DR method, at param wsir use a list with loadings in slot 2 (e.g PCA loadings) of dimension p * d
#' @param newdata matrix of new gene expression data to project into low-dimensional space
#'
#' @return matrix of low-dimensional representation of newdata gene expression data
#'
#' @examples
#' wsir_obj = wSIR(exprs = sample1_exprs, coords = sample1_coords)
#' sample2_low_dim_exprs = project_wSIR(wsir = wsir_obj, newdata = sample2_exprs)
#'
#' @export
project_wSIR = function(wsir, newdata) {

  ## Projecting new data
  # wsir is list object output from `wSIR`
  # newdata is a cells x features matrix of expression

  proj = newdata %*% wsir[[2]]
  return(proj)
}
