#' wSIR
#'
#' @description
#' A function to perform supervised dimension reduction of a gene expression matrix using coordinates dataframe as
#' response value. Incorporates weighting mechanism into SIR algorithm to generate low-dimensional representation of
#' the data and allow for projection of new single-cell gene expression data into low-dimensional space.
#'
#' @param X matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param samples sample ID of each cell. In total, must have length equal to the number of cells. For example, if
#' your dataset has 10000 cells, the first 5000 from sample 1 and the remaining 5000 from sample 2, you would write
#' samples = c(rep(1, 5000), rep(2, 5000)) to specify that the first 5000 cells are sample 1 and the remaining are sample 2.
#' Default is that all cells are from sample 1. Sample IDs can be of any format: for the previous example, you could write
#' samples = c(rep("sample 1", 5000), rep("sample 2", 5000)), and the result would be the same.
#' @param slices integer for number of slices on each axis of tissue. For example, slices = 4 creates 4 slices along
#' each spatial axis, yielding 4^2 = 16 tiles. Default is 8, suggested minimum of 3. Suggest to tune this parameter.
#' @param alpha integer to specify desired strength of spatial correlation. Suggest to tune this parameter on testing
#' dataset among e.g values c(0,2,4,8). alpha = 0 gives SIR implementation. Larger values give stronger spatial correlations.
#' @param maxDirections integer for upper limit on number of directions to include in low-dimensional space. Use if you
#' need less than a certain number of dimensions for downstream analyes
#' @param varThreshold numeric proportion of eigenvalues of variance in t(X_H) %*% W %*% X_H to retain. Default is 0.95.
#' Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#' @param optim_params logical, if you would like wSIR to automatically optimise parameters slices and alpha based on
#' either distance correlation or correlation of distances as evaluation metric. Default is TRUE. If your downstream
#' task is quite different to either of those metrics, then we suggest you tune those two parameters yourself using
#' your specific task and evaluation metric. In that case, determine your optimal slices and alpha values and then use
#' them in the relevant function, then setting optim_params = FALSE.
#' @param alpha_vals If you have optim_params = TRUE, then this is the values of alpha to optimise over in WSIR. 0
#' gives Sliced Inverse Regression (SIR) implementation, and larger values represent stronger spatial correlation.
#' Suggest to use integers for interpretability, but can use non-integers. Values must be non-negative.
#' @param slice_vals If you have optim_params = TRUE, then this is the values of slices to optimise over in WSIR.
#' Suggest maximum value in the vector to be no more than around \eqn{\sqrt{n/20}}, as this upper bound ensures an
#' average of at least 10 cells per tile in the training set.
#' @param metric If optim_params = TRUE, this is the evaluation metric to use for parameter tuning. String,
#' either "DC" to use distance correlation or "CD" to use correlation of distances. Default is "DC".
#' @param nrep If optim_params = TRUE, this is the integer for the number of train/test splits of the data to
#' perform during optimisation of parameters slices and alpha.
#'
#' @return wSIR object which is a list containing 5 (named) positions. 1) scores matrix containing low-dimensional
#' representation of X from wSIR algorithm. Dimension n * d, where d is chosen to include at least varThreshold proportion of variance.
#' 2) directions matrix containing WSIR directions, dimension p * d. Use this to project new data into low-dimensional space via X_new %*% directions.
#' 3) estd integer stating how many dimensions in the computed low-dimensional space are needed to account for varThreshold
#' proportion of variance. Same as number of dimensions in scores.
#' 4) W matrix weight matrix from cells_weight_matrix2 function
#' 5) evalues vector containing p eigenvalues of t(X_H) %*% W %*% X_H. varThreshold parameter works on these evalues,
#' such that e.g the first j directions are included if the sum of the first j evalues equals 0.95% of the sum of all evalues.
#'
#' @examples
#' data(MouseData)
#' wsir_obj = wSIR(X = sample1_exprs,
#'   coords = sample1_coords,
#'   optim_params = TRUE,
#'   alpha_vals = c(0,2,4),
#'   slice_vals = c(3,6,10)) # create wsir object
#'
#' @export
wSIR = function(X,
                coords,
                samples = rep(1, nrow(coords)),
                slices = 8,
                alpha = 4,
                maxDirections = 50,
                varThreshold = 0.95,
                optim_params = TRUE,
                alpha_vals = c(0,1,2,4,8,12),
                slice_vals = c(3,5,7,10,15,20),
                metric = "DC",
                nrep = 5) {

  if (optim_params) {
    optim_obj = exploreWSIRParams(exprs = X,
                                  coords = coords,
                                  samples = samples,
                                  alpha_vals = alpha_vals,
                                  slice_vals = slice_vals,
                                  varThreshold = varThreshold,
                                  maxDirections = maxDirections,
                                  metric = metric,
                                  nrep = nrep)
    alpha = optim_obj$best_alpha
    slices = optim_obj$best_slices
  }

  wsir_obj <- wSIRSpecifiedParams(X = X,
                                  coords = coords,
                                  samples = samples,
                                  alpha = alpha,
                                  slices = slices,
                                  maxDirections = maxDirections,
                                  varThreshold = varThreshold)
  return(wsir_obj)
}
