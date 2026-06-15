context("Test projectWSIR")

library(wSIR)
n <- 100
p <- 10
slice_vals <- c(5,8,10)
alpha_vals <- c(0,4,8)
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

wsir_out <- wSIR(X = x,
                 coords = coords)

new_n <- 80
new_x <- matrix(rnorm(n = new_n*p), nrow = new_n, ncol = p)

newdata_low_dim_exprs <- projectWSIR(wsir = wsir_out, 
                                     newdata = new_x)

expect_equal(nrow(newdata_low_dim_exprs), new_n)
expect_equal(ncol(newdata_low_dim_exprs), wsir_out$estd)