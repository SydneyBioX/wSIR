#' wSIROptimisation
#'
#' @description
#' This function is used to calculate the validation metric on which the best
#' tuning parameters are selected
#'
#' @param exprs_train matrix of normalised gene expression on the training data.
#' @param coords_train dataframe containing spatial positions of cells in 2D
#' space on the training data. Dimension n * 2. Column names must be
#' c("x", "y").
#' @param exprs_test matrix of normalised gene expression on the validation
#' data.
#' @param coords_train dataframe containing spatial positions of cells in 2D
#' space on the  validation data. Dimension n * 2. Column names must be
#' c("x", "y").
#' @param samples_train sample ID of each cell in the training data. In total,
#' must have length equal to the number of cells. For example, if
#' your dataset has 10000 cells, the first 5000 from sample 1 and the remaining
#' 5000 from sample 2, you would write
#' samples = c(rep(1, 5000), rep(2, 5000)) to specify that the first 5000 cells
#' are sample 1 and the remaining are sample 2.
#' Default is that all cells are from sample 1. Sample IDs can be of any
#' format: for the previous example, you could write
#' samples = c(rep("sample 1", 5000), rep("sample 2", 5000)), and the result
#' would be the same.
#' @param slices integer for number of slices on each axis of tissue. For
#' example, slices = 4 creates 4 slices along each spatial axis, yielding
#' 4^2 = 16 tiles.
#' Default is 8, suggested minimum of 3. Suggest to tune this parameter
#' using exploreWSIRParams() function.
#' @param alpha integer to specify desired strength of spatial correlation.
#' alpha = 0 gives SIR implementation. Larger values give stronger spatial
#' correlations.
#' Default value is 4, at which weight matrix value for a pair of tiles is
#' inversely proportional to the physical distance between them. Suggest
#' to tune this
#' parameter using exploreWSIRParams() function.
#' @param evalmetrics evaluation metrics to use for parameter tuning.
#' String, options are any or all of: "DC" to use distance
#' correlation; "CD" to use correlation of distances; "ncol" to use number
#' of columns in low-dimensional embedding. Default is all three,
#' specified by metrics = c("DC", "CD", "ncol").
#' @param ... arguments passed to wSIRSpecifiedParams
#'
#' @return Average metric value for the selected metric(s) over each
#' train/test split.
#'
#' @importFrom stats cor
#' @importFrom distances distances
#' @importFrom Rfast dcor
#'
#' @keywords internal

wSIROptimisation <- function(exprs_train,
                             coords_train,
                             exprs_test,
                             coords_test,
                             samples_train,
                             slices,
                             alpha,
                             evalmetrics = c("CD","DC","ncol"),
                             ...) {

    results <- NULL
    wsir_obj <- wSIRSpecifiedParams(X = exprs_train,
                                    coords = coords_train,
                                    samples = samples_train,
                                    slices = slices,
                                    alpha = alpha,
                                    ...)
   projected_test <- projectWSIR(wsir = wsir_obj, newdata = exprs_test)

    if ("CD" %in% evalmetrics) {
        # Replace dist by distances
        d1 <- as.matrix(distances::distances(projected_test))
        d2 <- as.matrix(distances::distances(coords_test))
        k1 <- .subsetLowerTri(as.matrix(d1))
        k2 <- .subsetLowerTri(as.matrix(d2))
        current_cd <- .spearman_correlation(k1, k2)
        results <- c(results, cd = current_cd)
    }
    if ("DC" %in% evalmetrics) {
        current_dc <- Rfast::bcdcor(as.matrix(projected_test),
                                    as.matrix(coords_test))
        results <- c(results, dc = current_dc)

    }
    if ("ncol" %in% evalmetrics) {
      current_ncol <- ncol(projected_test)
      ncol_vals <- current_ncol
      results <- c(results, ncol = ncol_vals)
    }

    return(results)
}
