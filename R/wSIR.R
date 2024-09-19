#' wSIR
#'
#' @description
#' A function to perform supervised dimension reduction of a gene expression
#' matrix using coordinates dataframe as
#' response value. Incorporates weighting mechanism into SIR algorithm to
#' generate low-dimensional representation of
#' the data and allow for projection of new single-cell gene expression data
#' into low-dimensional space.
#'
#' @param X matrix containing normalised gene expression data including n
#' cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D
#' space. Dimension n * 2. Column names must be c("x", "y").
#' @param optim_params logical default FALSE. If you would like wSIR to
#' automatically optimise parameters slices and alpha based on
#' either distance correlation or correlation of distances as evaluation
#' metric. If your downstream
#' task is quite different to either of those metrics, then we suggest you
#' tune those two parameters yourself using
#' your specific task and evaluation metric. In that case, determine your
#' optimal slices and alpha values and then use
#' them in the relevant function, then setting optim_params = FALSE.
#' @param alpha_vals If you have optim_params = TRUE, then this is the
#' values of alpha to optimise over in wSIR. 0
#' gives Sliced Inverse Regression (SIR) implementation, and larger values
#' represent stronger spatial correlation.
#' Suggest to use integers for interpretability, but can use non-integers.
#' Values must be non-negative.
#' @param slice_vals If you have optim_params = TRUE, then this is the values
#' of slices to optimise over in wSIR.
#' Suggest maximum value in the vector to be no more than around
#' \eqn{\sqrt{n/20}}, as this upper bound ensures an
#' average of at least 10 cells per tile in the training set.
#' @param metric If optim_params = TRUE, this is the evaluation metric to use
#' for parameter tuning. String,
#' either "DC" to use distance correlation or "CD" to use correlation of
#' distances. Default is "DC".
#' @param nrep If optim_params = TRUE, this is the integer for the number of
#' train/test splits of the data to
#' perform during optimisation of parameters slices and alpha.
#' @param verbose logical (default FALSE) whether progress messages should be
#' printed
#' @param ... arguments passing to `exploreWSIRParams` and `wSIRSpecifiedParams`
#'
#' @return wSIR object which is a list containing 5 (named) positions.
#' 1) scores matrix containing low-dimensional
#' representation of X from wSIR algorithm. Dimension `n * d`, where d is
#' chosen to include at least varThreshold proportion of variance.
#' 2) directions matrix containing wSIR directions, dimension `p * d`. Use
#' this to project new data into low-dimensional space via X_new %*% directions.
#' 3) estd integer stating how many dimensions in the computed
#' low-dimensional space are needed to account for varThreshold
#' proportion of variance. Same as number of dimensions in scores.
#' 4) W matrix weight matrix from cells_weight_matrix2 function
#' 5) evalues vector containing p eigenvalues of `t(X_H) %*% W %*% X_H`.
#' varThreshold parameter works on these evalues,
#' such that e.g the first j directions are included if the sum of the
#' first j evalues equals 0.95% of the sum of all evalues.
#'
#' @useDynLib wSIR
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' data(MouseData)
#' wsir_obj = wSIR(X = sample1_exprs,
#'   coords = sample1_coords,
#'   optim_params = TRUE,
#'   alpha_vals = c(0,2,4),
#'   slice_vals = c(3,6,10),
#'   metric = "DC",
#'   nrep = 1) # create wsir object
#'
#' @export
wSIR <- function(X,
                 coords,
                 optim_params = FALSE,
                 alpha_vals = c(0,1,2,4,8,12),
                 slice_vals = c(3,5,7,10,15,20),
                 metric = "DC",
                 nrep = 5,
                 verbose = FALSE,
                 ...) {

  if (optim_params) {
    if (verbose) message("Optimising parameters...")
    optim_obj <- exploreWSIRParams(X = X,
                                   coords = coords,
                                   alpha_vals = alpha_vals,
                                   slice_vals = slice_vals,
                                   metric = metric,
                                   nrep = nrep,
                                   ...)
    alpha <- optim_obj$best_alpha
    slices <- optim_obj$best_slices
    if (verbose) message("Optimising parameters... complete!")
  }

  if (verbose) message("Fitting wSIR model...")
  wsir_obj <- wSIRSpecifiedParams(X = X,
                                  coords = coords,
                                  ...
  )
  if (verbose) message("Fitting wSIR model... complete!")

  return(wsir_obj)
}
