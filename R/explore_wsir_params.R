#' explore_wsir_params function
#'
#' @description
#' This function is used to select the optimal values for parameters slices and alpha in weighted sliced inverse regression
#' based on your provided gene expression data and corresponding spatial coordinates. For a given evaluation metric,
#' it will visualise the performance of WSIR with varying function values based on your data.
#' NOTE: as of 4/7/24, this function won't work due to dependency on functions not in this repo. To be fixed sometime.
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
#' @param metric evaluation metric to use for parameterr tuning. String, either "DC" to use distance correlation or "CD" to use
#' correlation of distances. Default is "DC".
#' @param nrep integer for the number of train/test splits of the data to perform.
#' @param print_progress logical If TRUE, prints the current values of slices and alpha as the tuning gets up to performing WSIR
#' with each value. If FALSE, then no progress updates.
#'
#' @return Plot showing the average metric value across the nrep iterations for every combination of parameters slices and alpha.
#' Larger circles for a slices/alpha combination indicates better performance for that pair of values. Suggest to find the pair with
#' the highest metric value and use those parameter values for your task.
#'
#' @examples
#' plot = explore_wsir_params(exprs = sample1_exprs,
#' coords = sample1_coords,
#' alpha_vals = c(0,2,4,8),
#' slice_vals = c(3,6,10))
#' plot
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
                               print_progress = TRUE) {

  if (metric == "DC") {
    distCor = TRUE
    corDist = FALSE
  } else if (metric == "CD") {
    corDist = TRUE
    distCor = FALSE
  }

  metric_vals <- c()

  for (alpha in alpha_vals) {
    if (print_progress) {
      print(paste("current alpha:", alpha))
    }
    for (slices in slice_vals) {
      if (print_progress) {
        print(paste("current slices:", slices))
      }
      metric_current <- analysis_all(WSIR = TRUE,
                                   exprs = exprs,
                                   coords = coords,
                                   slices = slices,
                                   alpha = alpha,
                                   distCor = distCor,
                                   corDist = corDist,
                                   directions = maxDirections,
                                   varThreshold = varThreshold,
                                   nrep = nrep,
                                   plot_show = FALSE)$dataframe[,2] %>% mean() # take average metric value over nrep iterations
      metric_vals <- metric_vals %>% append(metric_current)
    }
  }

  res_df <- matrix(NA, nrow = length(alpha_vals)*length(slice_vals), ncol = 3) %>% as.data.frame()
  colnames(res_df) <- c("alpha", "slices", "metric")

  res_df$alpha <- vec_rep_each(alpha_vals, length(slice_vals)) %>% as.factor()
  res_df$slices <- rep(slice_vals, length(alpha_vals)) %>% as.factor()
  res_df$metric <- metric_vals

  plot <- ggplot(data = res_df, aes(x = alpha, y = slices, size = metric)) +
    geom_point() +
    theme_classic() +
    ggtitle("Metric value for different parameter combinations (nrep iterations of train/test split)")

  return(plot)

}
