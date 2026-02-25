context("Test findTopGenes")

library(wSIR)

# check that the dimension of findTopGenes$genes dataframe is correct dimension
n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

out <- wSIR(X = x,
            coords = coords)
nhighest <- 3
ndirs <- 4
find_out <- wSIR::findTopGenes(WSIR = out,
                               highest = nhighest,
                               dirs = ndirs)
dim_find_out <- nrow(find_out$genes)
expect_equal(dim_find_out, nhighest*ndirs)
