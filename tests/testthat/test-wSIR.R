context("Test wSIR")

library(wSIR)
# This is for testing for an error with invalid parameters

n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

slices_good <- 3
slices_bad <- 3.5
alpha_good <- 4
alpha_bad <- -4
maxd_good <- 30
maxd_bad <- 30.5
varThreshold_good <- 0.95
varThreshold_bad <- 1.01

# slices
expect_error(out <- wSIR(X = x,
                         coords = coords,
                         slices = slices_bad,
                         alpha = alpha_good,
                         maxDirections = maxd_good,
                         varThreshold = varThreshold_good))
# alpha
expect_error(out <- wSIR(X = x,
                         coords = coords,
                         slices = slices_good,
                         alpha = alpha_bad,
                         maxDirections = maxd_good,
                         varThreshold = varThreshold_good))
# maxDirections
expect_error(out <- wSIR(X = x,
                         coords = coords,
                         slices = slices_good,
                         alpha = alpha_good,
                         maxDirections = maxd_bad,
                         varThreshold = varThreshold_good))
# varThreshold
expect_error(out <- wSIR(X = x,
                         coords = coords,
                         slices = slices_good,
                         alpha = alpha_good,
                         maxDirections = maxd_good,
                         varThreshold = varThreshold_bad))

## New test: test to expect an error if nrow(X), nrow(coords), length(samples) not all equal
# Setup: expect error if any or all of nrow(X), nrow(coords), length(samples) not same
n1 <- 100
n2 <- 101
n3 <- 102
p <- 10
x1 <- matrix(rnorm(n = n1*p), nrow = n1, ncol = p)
x2 <- matrix(rnorm(n = n2*p), nrow = n2, ncol = p)
x3 <- matrix(rnorm(n = n3*p), nrow = n3, ncol = p)
coords1 <- matrix(runif(n <- n1*2), nrow = n1, ncol = 2)
coords2 <- matrix(runif(n <- n2*2), nrow = n2, ncol = 2)
coords3 <- matrix(runif(n <- n3*2), nrow = n3, ncol = 2)
samples1 <- sample(c(1,2), size = n1, replace = TRUE)
samples2 <- sample(c(1,2), size = n2, replace = TRUE)
samples3 <- sample(c(1,2), size = n3, replace = TRUE)

# length(samples) unequal
expect_error(out <- wSIR(X = x1,
                         coords = coords1,
                         samples = samples2))
# nrow(coords) unequal
expect_error(out <- wSIR(X = x1,
                         coords = coords2,
                         samples = samples1))
# nrow(X) unequal
expect_error(out <- wSIR(X = x2,
                         coords = coords1,
                         samples = samples1))
# all unequal
expect_error(out <- wSIR(X = x1,
                         coords = coords2,
                         samples = samples3))










