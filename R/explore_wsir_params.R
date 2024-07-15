#' explore_wsir_params function
#'
#' @description
#' This function is used to select the optimal values for parameters slices and alpha in weighted sliced inverse regression
#' based on your provided gene expression data and corresponding spatial coordinates. For a given evaluation metric,
#' it will visualise the performance of WSIR with varying function parameters based on your data, and return the optimal pair.
#' This pair of slices and alpha can be used for your downstream tasks.
#'
#' @param exprs matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param alpha_vals vector of numbers as the values of parameter alpha to use in WSIR. 0 gives Sliced Inverse Regression
#' (SIR) implementation, and larger values represent stronger spatial correlation. Suggest to use integers for interpretability,
#' but can use non-integers. Values must be non-negative.
#' @param slice_vals vector of integers as the values of parameter slices to use in WSIR. Suggest maximum value in the vector to
#' be no more than around \eqn{\sqrt{n/20}}, as this upper bound ensures an average of at least 10 cells per tile in the training set.
#' @param varThreshold numeric proportion of variance in \code{t(X_H) \%*\% W \%*\% X_H} to retain. Must be between 0 and 1. Default is 0.95.
#' Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#' @param maxDirections integer for the maximum number of directions to include in the low-dimenensional embedding. Default is 50.
#' @param metric evaluation metric to use for parameter tuning. String, either "DC" to use distance correlation or "CD" to use
#' correlation of distances. Default is "DC".
#' @param nrep integer for the number of train/test splits of the data to perform.
#' @param verbose logical If TRUE, prints the current values of slices and alpha as the tuning gets up to performing WSIR
#' with each value. If FALSE, then no progress updates. Default is FALSE.
#'
#' @return List with four slots, named "plot" and "message".
#' 1) "Plot" shows the average metric value across the nrep iterations for every combination of parameters slices and alpha.
#' Larger circles for a slices/alpha combination indicates better performance for that pair of values.
#' 2) "message" tells you the parameter combination with highest metric value.
#' 3) "best_alpha" returns the integer for the best alpha values among the values that were tested.
#' 4) "best_slices" returns the integer for the best slices value among the values that were tested.
#'
#' @examples
#' data(MouseData)
#' explore_params = explore_wsir_params(exprs = sample1_exprs, coords = sample1_coords, alpha_vals = c(0,2,4,8), slice_vals = c(3,6,10))
#' explore_params$plot
#' explore_params$message
#' wsir_obj = wSIR(X = sample1_exprs, coords = sample1_coords, alpha = best_alpha, slices = best_slices)
#'
#' @export
explore_wsir_params = function(exprs,
                               coords,
                               alpha_vals = c(0,1,2,4,8,12),
                               slice_vals = c(3,5,7,10,15,20),
                               varThreshold = 0.95,
                               maxDirections = 50,
                               metric = "DC",
                               nrep = 5,
                               verbose = FALSE) {

  metric_vals <- c()

  for (alpha in alpha_vals) {
    if (verbose) {
      print(paste("current alpha:", alpha))
    }
    for (slices in slice_vals) {
      if (verbose) {
        print(paste("current slices:", slices))
      }
      metric_current <- metric_eval_wsir_optim(exprs = exprs,
                                               coords = coords,
                                               alpha = alpha,
                                               slices = slices,
                                               varThreshold = varThreshold,
                                               maxDirections = maxDirections,
                                               metric = metric,
                                               nrep = nrep)
      metric_vals <- metric_vals %>% append(metric_current)
    }
  }

  res_df <- matrix(NA, nrow = length(alpha_vals)*length(slice_vals), ncol = 3) %>% as.data.frame()
  colnames(res_df) <- c("alpha", "slices", "metric")

  res_df$alpha <- vec_rep_each(alpha_vals, length(slice_vals))
  res_df$slices <- rep(slice_vals, length(alpha_vals))
  res_df$metric <- metric_vals

  best_alpha = res_df$alpha[which.max(res_df$metric)]
  best_slices = res_df$slices[which.max(res_df$metric)]

  res_df$alpha <- res_df$alpha %>% as.factor()
  res_df$slices <- res_df$slices %>% as.factor()

  message = paste0("Optimal (alpha, slices) pair: (", best_alpha, ", ", best_slices, ")")

  plot <- ggplot(data = res_df, aes(x = alpha, y = slices, size = metric)) +
    geom_point() +
    theme_classic() +
    ggtitle(paste0("Metric value for different parameter combinations (",nrep, " iterations of train/test split)"))

  return(list(plot = plot,
              message = message,
              best_alpha = best_alpha,
              best_slices = best_slices))

}
