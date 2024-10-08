#' exploreWSIRParams function
#'
#' @description
#' This function is used to select the optimal values for parameters slices
#' and alpha in weighted sliced inverse regression
#' based on your provided gene expression data and corresponding spatial
#' coordinates. For a given evaluation metric,
#' it will visualise the performance of WSIR with varying function parameters
#' based on your data, and return the optimal pair.
#' This pair of slices and alpha can be used for your downstream tasks.
#'
#' @param X matrix containing normalised gene expression data including n
#' cells and p genes, dimension n * p.
#' @param coords dataframe containing spatial positions of n cells in 2D space.
#' Dimension n * 2. Column names must be c("x", "y").
#' @param samples sample ID of each cell. In total, must have length equal to
#' the number of cells. For example, if
#' your dataset has 10000 cells, the first 5000 from sample 1 and the remaining
#' 5000 from sample 2, you would write
#' samples = c(rep(1, 5000), rep(2, 5000)) to specify that the first 5000 cells
#' are sample 1 and the remaining are sample 2.
#' Default is that all cells are from sample 1. Sample IDs can be of any
#' format: for the previous example, you could write
#' samples = c(rep("sample 1", 5000), rep("sample 2", 5000)), and the result
#' would be the same.
#' @param alpha_vals vector of numbers as the values of parameter alpha to use
#' in WSIR. 0 gives Sliced Inverse Regression
#' (SIR) implementation, and larger values represent stronger spatial
#' correlation. Suggest to use integers for interpretability,
#' but can use non-integers. Values must be non-negative.
#' @param slice_vals vector of integers as the values of parameter slices to
#' use in WSIR. Suggest maximum value in the vector to
#' be no more than around \eqn{\sqrt{n/20}}, as this upper bound ensures an
#' average of at least 10 cells per tile in the training set.
#' @param metric character single value. Evaluation metric to use for parameter
#' tuning to select optimal parameter combination. String, use "DC" to
#' use distance correlation, "CD" to use correlation of distances, or "ncol"
#' for the number of dimensions in the low-dimensional
#' embedding. Default is "DC".
#' @param nrep integer for the number of train/test splits of the data to
#' perform.
#' @param param parallel computing setup for bplapply from BiocParallel
#' package. Default is to use a single core, hence
#' default value is MulticoreParam(workers = 1)
#' @param ... arguments passed on to wSIROptimisation
#'
#' @return List with five slots, named "plot", "message", "best_alpha",
#' "best_slices" and "results_dataframe".
#' 1) "plot" shows the average metric value across the nrep iterations for
#' every combination of parameters slices and alpha.
#' Larger circles for a slices/alpha combination indicates better performance
#' for that pair of values. There is one panel per
#' evaluation metric selected in "metrics" argument.
#' 2) "message" tells you the parameter combination with highest metric value
#' according to selected metric.
#' 3) "best_alpha" returns the integer for the best alpha values among the
#' values that were tested according to selected metric.
#' 4) "best_slices" returns the integer for the best slices value among the
#' values that were tested according to selected metric.
#' 5) "results_dataframe" returns the results dataframe used to create "plot".
#' This dataframe has length(alpha_vals)*length(slice_vals) rows,
#' where one is for each combination of parameters slices and alpha. There is
#' one column for "alpha", one for "slices" and one
#' for each of the evaluation metrics selected in "metrics" argument. Column
#' "alpha" includes the value for parameter alpha,
#' column "slices" includes the value for parameter slices, and each metric
#' column includes the value for the specified metric,
#' which is either Distance Correlation ("DC"), Correlation of Distances
#' ("CD"), or number of columns in low-dimensional embedding ("ncol").
#'
#' @examples
#' data(MouseData)
#' explore_params = exploreWSIRParams(X = sample1_exprs,
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
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle
#' @importFrom vctrs vec_rep_each
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom stringr word
#'
#' @export
exploreWSIRParams <- function(X,
                              coords,
                              samples = rep(1, nrow(coords)),
                              alpha_vals = c(0,2,4,10),
                              slice_vals = c(5,10,15,20),
                              metric = "DC",
                              nrep = 5,
                              param = MulticoreParam(workers = 1),
                              ...
) {

  # vector of all parameter combinations
  param_combinations <- expand.grid(slice = slice_vals, alpha = alpha_vals,
                                    metric = NA)


  #  # perform bplapply over that list of (pairs of) combinations
  #  metric_vals_list = bplapply(param_combinations, function(parameter_pair) {
  #    current_slices = as.numeric(word(parameter_pair, 1, sep = ","))
  #    current_alpha = as.numeric(word(parameter_pair, 2, sep = ","))
  #    optim_result = wSIROptimisation(exprs = exprs,
  #                                    coords = coords,
  #                                    alpha = current_alpha,
  #                                    slices = current_slices,
  #                                   varThreshold = varThreshold,
  #                                   maxDirections = maxDirections,
  #                                   metrics = metrics,
  #                                   nrep = nrep)
  #   return(optim_result)
  # },
  # BPPARAM = param)

  # Create pre-specified random splits of data, each columns
  # corresponding to one split
  index_rep <- matrix(
    sample(c(TRUE, FALSE), nrow(X)*nrep, replace = TRUE),
    nrow = nrow(X), ncol = nrep
  )
  # create training and test set from each column index
  split_list <- apply(index_rep, 2,  function(keep){
    X_train <- X[keep,]
    coords_train <- coords[keep,]
    samples_train <- samples[keep]
    X_test <- X[!keep,]
    coords_test <- coords[!keep,]
    list(X_train, coords_train, X_test,
         coords_test,  samples_train )
  })
  nElements <- 5
  result <- lapply(seq_len(nElements),
                   function(i) lapply(split_list, "[[", i))

  for (ii in seq_len(nrow(param_combinations))) {
    cv_scores <- mapply(function(X_train, coords_train, X_test,
                                 coords_test,  samples_train){
      wSIROptimisation(as.matrix(X_train),
                       coords_train, as.matrix(X_test),
                       coords_test,
                       samples_train, param_combinations$slice[ii],
                       param_combinations$alpha[ii],
                       evalmetrics = metric,
                       ...)
    }, result[[1]], result[[2]], result[[3]], result[[4]], result[[5]])

    param_combinations$metric[ii] <- mean(cv_scores, na.rm = TRUE)
  }

  res_df <- param_combinations
  best_metric_index <- which.max(res_df[, "metric"])
  best_alpha <- res_df[best_metric_index, "alpha"]
  best_slices <- res_df[best_metric_index, "slice"]

  message <- paste0("Optimal (alpha, slices) pair: (",
                    best_alpha, ", ", best_slices, ")")

  plot <- ggplot2::ggplot(res_df, mapping = aes(
    x = .data$alpha, y = .data$slice, size = .data$metric)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(
      paste0("Metric value for different parameter combinations (",
             nrep, " iterations of train/test split)"))

  return(list(plot = plot,
              message = message,
              best_alpha = best_alpha,
              best_slices = best_slices,
              results_dataframe = res_df))

}
