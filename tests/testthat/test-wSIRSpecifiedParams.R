context("Test wSIRSpecifiedParams")

library(wSIR)
n <- 100
p <- 10
slices <- 3
alpha <- 4
maxd <- 30
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

# test for maximum directions being followed
out <- wSIR:::wSIRSpecifiedParams(X = x,
                                  coords = coords,
                                  slices = slices,
                                  alpha = alpha,
                                  maxDirections = 30)
dim_scores <- ncol(out$scores)
dim_dirs <- nrow(out$directions)
expect_lte(dim_scores, maxd)
expect_lte(dim_dirs, maxd)