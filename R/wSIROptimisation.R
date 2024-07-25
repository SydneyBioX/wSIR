#' wSIROptimisation
#'
#' @description
#' This function is used for parameter optimisation in WSIR. with provided parameter values, exprs and coords, this function splits the data into equal train/test
#' splits a given number of times (nrep times). It computes the selected evaluation metric (distance correlation or correlation of distances) on the current test
#' set and returns the average value for that metric over all train/test splits.
#'
#' @param exprs matrix of normalised gene expression data for n genes across p cells.
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param samples sample ID of each cell. In total, must have length equal to the number of cells. For example, if
#' your dataset has 10000 cells, the first 5000 from sample 1 and the remaining 5000 from sample 2, you would write
#' samples = c(rep(1, 5000), rep(2, 5000)) to specify that the first 5000 cells are sample 1 and the remaining are sample 2.
#' Default is that all cells are from sample 1. Sample IDs can be of any format: for the previous example, you could write
#' samples = c(rep("sample 1", 5000), rep("sample 2", 5000)), and the result would be the same.
#' @param slices integer for number of slices on each axis of tissue. For example, slices = 4 creates 4 slices along each spatial axis, yielding 4^2 = 16 tiles.
#' Default is 8, suggested minimum of 3. Suggest to tune this parameter using exploreWSIRParams() function.
#' @param alpha integer to specify desired strength of spatial correlation. alpha = 0 gives SIR implementation. Larger values give stronger spatial correlations.
#' Default value is 4, at which weight matrix value for a pair of tiles is inversely proportional to the physical distance between them. Suggest to tune this
#' parameter using exploreWSIRParams() function.
#' @param maxDirections integer for the maximum number of directions to include in the low-dimenensional embedding. Default is 50.
#' @param varThreshold numeric proportion of variance in \code{t(X_H) \%*\% W \%*\% X_H} to retain. Must be between 0 and 1. Default is 0.95.
#' Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#' @param metric evaluation metric to use for parameterr tuning. String, either "DC" to use distance correlation or "CD" to use
#' correlation of distances. Default is "DC".
#' @param nrep integer for the number of train/test splits of the data to perform.
#'
#' @return Average metric value for the given metric over each train/test split.
#'
#' @importFrom stats cor
#' @importFrom stats dist
#'
#' @keywords internal

wSIROptimisation = function(exprs,
                            coords,
                            samples = rep(1, nrow(coords)),
                            slices,
                            alpha,
                            maxDirections,
                            varThreshold,
                            metric,
                            nrep = 3) {
  metric_vals <- rep(0, nrep)
  for (i in 1:nrep) {
    keep <- sample(c(TRUE, FALSE), nrow(exprs), replace = TRUE)
    exprs_train <- exprs[keep,]
    coords_train <- coords[keep,]
    samples_train <- samples[keep]

    exprs_test <- exprs[!keep,]
    coords_test <- coords[!keep,]

    wsir_obj = wSIRSpecifiedParams(X = exprs_train,
                                   coords = coords_train,
                                   samples = samples_train,
                                   slices = slices,
                                   alpha = alpha,
                                   maxDirections = maxDirections,
                                   varThreshold = varThreshold)
    projected_test = projectWSIR(wsir = wsir_obj, newdata = exprs_test)

    if (metric == "CD") {
      current_metric = cor(subsetLowerTri(dist(projected_test)),
                           subsetLowerTri(dist(coords_test)),
                           method = "spearman",
                           use = "pairwise.complete")
    } else if (metric == "DC") {
      current_metric = Rfast::dcor(x = as.matrix(projected_test),
                                   y = coords_test)$dcor
    }
    metric_vals[i] = current_metric
  }
  avg_metric = mean(metric_vals)
  return(avg_metric)
}
