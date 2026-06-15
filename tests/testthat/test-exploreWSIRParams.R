context("Test exploreWSIRParams")

library(wSIR)
n <- 100
p <- 10
slice_vals <- c(5,8,10)
alpha_vals <- c(0,4,8)
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

# check that the returned best alpha and slices are actually best
explore_params = exploreWSIRParams(X = x,
                                   coords = coords,
                                   optim_alpha = alpha_vals,
                                   optim_slices = slice_vals)

out_df <- explore_params$results_dataframe
best_metric_index <- which.max(out_df[,"metric"])
best_slices_checked <- out_df$slices[best_metric_index]
best_alpha_checked <- out_df$alpha[best_metric_index]

expect_equal(best_slices_checked, explore_params$best_slices)
expect_equal(best_alpha_checked, explore_params$best_alpha)