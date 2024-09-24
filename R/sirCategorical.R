#' sirCategorical
#'
#' @description
#' This function performs WSIR based on provided gene expression data, tile
#' allocations, and weight matrix.
#'
#' @param X matrix of normalised gene expression data for n genes across p
#' cells.
#' @param Y dataframe with 1 column named "coordinate" that is the tile
#' allocation for each cell. There should
#' be up to slices^2 unique tile IDs in this column.
#' @param W Weight matrix created by createWeightMatrix. Entry (i,j)
#' represents the spatial correlation level
#' between tiles i and j. The diagonal values should be all 1. If not
#' provided, SIR implementation will be used.
#'
#' @return list of outputs with 5 named slots. They are the same as the
#' output of the wSIR function: this is
#' the final step in the wSIR function.
#'
#' @keywords internal

sirCategorical <- function(X,
                           Y,
                           # maxDirections = 50,
                           W = NULL,
                           ...
) {

  # do the transformation (QR scaling method)
  n <- nrow(X)
  X <- as.matrix(X)
  RandZ <- .computeRandZ(X)
  R <- RandZ[[1]]
  Z <- RandZ[[2]]

  sliced_data <- slicerCategorical(X = Z, Y = Y)

  if (is.null(W)) {
    W <- base::diag(table(Y$coordinate), ncol = nrow(sliced_data))/nrow(Y)
  }

  pc_dirs <- sirPCA(sliced_data,
                    W = W,
                    ...
  )

  betas <- base::backsolve(R, pc_dirs$evectors)
  betas <- as.matrix(betas)
  betas <- apply(betas, 2, function(x) x/sqrt(sum(x^2)))
  rownames(betas) <- colnames(X)
  final_XB <- .matMultArma(X, betas)

  return(list(scores = final_XB,
              directions = betas,
              estd = pc_dirs[[2]],
              W = W,
              evalues = pc_dirs$evalues))
}
