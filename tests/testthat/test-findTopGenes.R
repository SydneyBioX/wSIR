context("Test findTopGenes")

library(wSIR)

# check that input X matrix has column names (otherwise no genes to find)
n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

out <- wSIR(X = x,
            coords = coords)
nhighest <- 3
ndirs <- 4
expect_error(find_out <- wSIR::findTopGenes(WSIR = out,
                                            highest = nhighest,
                                            dirs = c(1:ndirs)))

# check that the dimension of findTopGenes$genes dataframe is correct dimension
# make sure to give x column names now
n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
colnames(x) <- paste0("gene", c(1:p))
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

out <- wSIR(X = x,
            coords = coords)
nhighest <- 3
ndirs <- 4
find_out <- wSIR::findTopGenes(WSIR = out,
                               highest = nhighest,
                               dirs = c(1:ndirs))
dim_find_out <- nrow(find_out$genes)
expect_equal(dim_find_out, nhighest*ndirs)
