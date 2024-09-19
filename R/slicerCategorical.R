#' slicerCategorical
#'
#' @description
#' This function averages all columns in the X matrix for each category in the
#' single-column Y dataframe. It
#' is used to find the average position within each tile to create the weight
#' matrix, as well as to create the
#' tile means that are used in the eigendecomposition step of SIR/wSIR.
#'
#' @param X matrix or dataframe containing the columns to average
#' @param Y dataframe with a single column named 'coordinate'. This column
#' contains the categorical variable
#' that is the tile ID for each cell.
#'
#' @return dataframe containing the averages for each column within each
#' category in Y. Dimension h * p,
#' where h is the number of categories in Y and p is the number of columns in X.
#'
#'
#' @keywords internal

slicerCategorical <- function(X, Y) {

  avg_X <- do.call(rbind,lapply(split.data.frame(X, Y$coordinate), colMeans))

  return(avg_X)

}
