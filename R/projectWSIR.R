#' projectWSIR
#'
#' @description
#' function to project new gene expression data into low-dimensional space
#'
#' @param wsir wsir object that is usually the output of wSIR function. If you
#' want to project new data into low-dim space following a
#' different DR method, at param wsir use a list with matrix of loadings in
#' slot 2 (e.g PCA loadings) of dimension p * d
#' @param new_data matrix of new gene expression data to project into
#' low-dimensional space. Must have the same p columns
#' as the columns in X argument used to generate wsir.
#'
#' @return matrix of low-dimensional representation of newdata gene
#' expression data
#'
#' @examples
#' data(MouseData)
#' wsir_obj <- wSIR(X = sample1_exprs,
#'     coords = sample1_coords,
#'     optim_params = FALSE,
#'     alpha = 4,
#'     slices = 6)
#' sample2_low_dim_exprs <- projectWSIR(wsir = wsir_obj, 
#'     new_data = sample2_exprs)
#'
#' @export

projectWSIR <- function(wsir, new_data) {
    new_data <- as.matrix(new_data)
    proj <- .matMultArma(new_data, wsir[[2]])
    return(proj)
}
