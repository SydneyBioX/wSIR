#' metric_eval_wsir_optim
#'
#' @description
#' This function is used for parameter optimisation in WSIR. with provided parameter values, exprs and coords, this function splits the data into equal train/test 
#' splits a given number of times (nrep times). It computes the selected evaluation metric (distance correlation or correlation of distances) on the current test
#' set and returns the average value for that metric over all train/test splits.
#'
#' @param exprs matrix of normalised gene expression data for n genes across p cells. 
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param slices integer for number of slices on each axis of tissue. For example, slices = 4 creates 4 slices along each spatial axis, yielding 4^2 = 16 tiles. 
#' Default is 8, suggested minimum of 3. Suggest to tune this parameter using explore_wsir_params() function.
#' @param alpha integer to specify desired strength of spatial correlation. alpha = 0 gives SIR implementation. Larger values give stronger spatial correlations. 
#' Default value is 4, at which weight matrix value for a pair of tiles is inversely proportional to the physical distance between them. Suggest to tune this 
#' parameter using explore_wsir_params() function.
#' @param maxDirections integer for the maximum number of directions to include in the low-dimenensional embedding. Default is 50.
#' @param varThreshold numeric proportion of variance in \code{t(X_H) \%*\% W \%*\% X_H} to retain. Must be between 0 and 1. Default is 0.95.
#' Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#' @param metric evaluation metric to use for parameterr tuning. String, either "DC" to use distance correlation or "CD" to use
#' correlation of distances. Default is "DC".
#' @param nrep integer for the number of train/test splits of the data to perform.
#'
#' @return Average metric value for the given metric over each train/test split. 
#'
#' @keywords internal

metric_eval_wsir_optim = function(exprs,
                                  coords,
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
    
    exprs_test <- exprs[!keep,]
    coords_test <- coords[!keep,]
    
    wsir_obj = WSIR(X = exprs_train,
                    coords = coords_train,
                    slices = slices,
                    alpha = alpha,
                    maxDirections = maxDirections,
                    varThreshold = varThreshold)
    projected_test = project_wSIR(wsir = wsir_obj, newdata = exprs_test)
    
    if (metric == "CD") {
      current_metric = cor(subset_lower.tri(dist(projected_test)),
                           subset_lower.tri(dist(coords_test)),
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
