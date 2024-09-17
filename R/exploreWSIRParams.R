#' exploreWSIRParams function
#'
#' @description
#' This function is used to select the optimal values for parameters slices and alpha in weighted sliced inverse regression
#' based on your provided gene expression data and corresponding spatial coordinates. For a given evaluation metric,
#' it will visualise the performance of WSIR with varying function parameters based on your data, and return the optimal pair.
#' This pair of slices and alpha can be used for your downstream tasks.
#'
#' @param exprs matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param samples sample ID of each cell. In total, must have length equal to the number of cells. For example, if
#' your dataset has 10000 cells, the first 5000 from sample 1 and the remaining 5000 from sample 2, you would write
#' samples = c(rep(1, 5000), rep(2, 5000)) to specify that the first 5000 cells are sample 1 and the remaining are sample 2.
#' Default is that all cells are from sample 1. Sample IDs can be of any format: for the previous example, you could write
#' samples = c(rep("sample 1", 5000), rep("sample 2", 5000)), and the result would be the same.
#' @param alpha_vals vector of numbers as the values of parameter alpha to use in WSIR. 0 gives Sliced Inverse Regression
#' (SIR) implementation, and larger values represent stronger spatial correlation. Suggest to use integers for interpretability,
#' but can use non-integers. Values must be non-negative.
#' @param slice_vals vector of integers as the values of parameter slices to use in WSIR. Suggest maximum value in the vector to
#' be no more than around \eqn{\sqrt{n/20}}, as this upper bound ensures an average of at least 10 cells per tile in the training set.
#' @param varThreshold numeric proportion of variance in \code{t(X_H) \%*\% W \%*\% X_H} to retain. Must be between 0 and 1. Default is 0.95.
#' Select higher threshold to include more dimensions, lower threshold to include less dimensions.
#' @param maxDirections integer for the maximum number of directions to include in the low-dimensional embedding. Default is 50.
#' @param metrics evaluation metrics to display in the plots. In a vector, include any or all of "DC" for distance correlation,
#' "CD" for correlation of distances, or "ncol" for number of columns in the low-dimensional embedding. Default is all three, specified
#' by metrics = c("DC", "CD", "ncol").
#' @param metric evaluation metric to use for parameter tuning to select optimal parameter combination. String, use "DC" to
#' use distance correlation, "CD" to use correlation of distances, or "ncol" for the number of dimensions in the low-dimensional
#' embedding. Default is "DC".
#' @param nrep integer for the number of train/test splits of the data to perform.
#' @param param parallel computing setup for bplapply from BiocParallel package. Default is to use a single core, hence
#' default value is MulticoreParam(workers = 1)
#'
#' @return List with five slots, named "plot", "message", "best_alpha", "best_slices" and "results_dataframe".
#' 1) "plot" shows the average metric value across the nrep iterations for every combination of parameters slices and alpha.
#' Larger circles for a slices/alpha combination indicates better performance for that pair of values. There is one panel per
#' evaluation metric selected in "metrics" argument.
#' 2) "message" tells you the parameter combination with highest metric value according to selected metric.
#' 3) "best_alpha" returns the integer for the best alpha values among the values that were tested according to selected metric.
#' 4) "best_slices" returns the integer for the best slices value among the values that were tested according to selected metric.
#' 5) "results_dataframe" returns the results dataframe used to create "plot". This dataframe has length(alpha_vals)*length(slice_vals) rows,
#' where one is for each combination of parameters slices and alpha. There is one column for "alpha", one for "slices" and one
#' for each of the evaluation metrics selected in "metrics" argument. Column "alpha" includes the value for parameter alpha,
#' column "slices" includes the value for parameter slices, and each metric column includes the value for the specified metric,
#' which is either Distance Correlation ("DC"), Correlation of Distances ("CD"), or number of columns in low-dimensional embedding ("ncol").
#'
#' @examples
#' data(MouseData)
#' explore_params = exploreWSIRParams(exprs = sample1_exprs,
#'   coords = sample1_coords,
#'   alpha_vals = c(0,2,4,8),
#'   slice_vals = c(3,6,10))
#' explore_params$plot
#' explore_params$message
#' best_alpha = explore_params$best_alpha
#' best_slices = explore_params$best_slices
#' wsir_obj = wSIR(X = sample1_exprs,
#'   coords = sample1_coords,
#'   optim_params = FALSE,
#'   alpha = best_alpha,
#'   slices = best_slices)
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 ggtitle
#' @importFrom vctrs vec_rep_each
#' @importFrom BiocParallel bplapply
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stringr word
#'
#' @export
exploreWSIRParams = function(exprs,
                             coords,
                             samples = rep(1, nrow(coords)),
                             alpha_vals = c(0,2,4,10),
                             slice_vals = c(5,10,15,20),
                             varThreshold = 0.95,
                             maxDirections = 50,
                             metrics = c("DC", "CD", "ncol"),
                             metric = "DC",
                             nrep = 5,
                             param = MulticoreParam(workers = 1)) {

  # vector of all parameter combinations
  param_combinations = as.vector(outer(slice_vals, alpha_vals, paste, sep = ","))

  # perform bplapply over that list of (pairs of) combinations
  metric_vals_list = bplapply(param_combinations, function(parameter_pair) {
    current_slices = as.numeric(word(parameter_pair, 1, sep = ","))
    current_alpha = as.numeric(word(parameter_pair, 2, sep = ","))
    optim_result = wSIROptimisation(exprs = exprs,
                                    coords = coords,
                                    alpha = current_alpha,
                                    slices = current_slices,
                                    varThreshold = varThreshold,
                                    maxDirections = maxDirections,
                                    metrics = metrics,
                                    nrep = nrep)
    return(optim_result)
  },
  BPPARAM = param)

  metric_vals <- unlist(metric_vals_list)

  res_df <- matrix(NA, nrow = length(alpha_vals)*length(slice_vals)*length(metrics), ncol = 4) %>% as.data.frame()
  colnames(res_df) <- c("alpha", "slices", "metric", "value")

  res_df$alpha <- vec_rep_each(alpha_vals, length(slice_vals)*length(metrics))
  res_df$slices <- rep(vec_rep_each(slice_vals, length(metrics)), length(alpha_vals))
  res_df$metric <- rep(metrics, length(slice_vals)*length(alpha_vals))
  res_df$value <- metric_vals

  chosen_metric_df <- res_df[(res_df$metric == metric),]

  best_alpha = chosen_metric_df$alpha[which.max(chosen_metric_df$value)]
  best_slices = chosen_metric_df$slices[which.max(chosen_metric_df$value)]

  res_df$alpha <- res_df$alpha %>% as.factor()
  res_df$slices <- res_df$slices %>% as.factor()

  message = paste0("Optimal (alpha, slices) pair: (", best_alpha, ", ", best_slices, ")")

  plot <- ggplot(data = res_df, aes(x = alpha, y = slices, size = value)) +
    geom_point() +
    theme_classic() +
    ggtitle(paste0("Metric value for different parameter combinations (",nrep, " iterations of train/test split)")) +
    facet_wrap(~metric, scales = "free")

  return(list(plot = plot,
              message = message,
              best_alpha = best_alpha,
              best_slices = best_slices,
              results_dataframe = res_df))

}
