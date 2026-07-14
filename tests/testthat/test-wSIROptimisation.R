#context("Test wSIROptimisation")

library(wSIR)
n1 <- 100
n2 <- 50
p <- 10
x1 <- matrix(rnorm(n = n1*p), nrow = n1, ncol = p)
x2 <- matrix(rnorm(n = n2*p), nrow = n2, ncol = p)
colnames(x1) <- paste0("g", c(1:p))
colnames(x2) <- paste0("g", c(1:p))
coords1 <- data.frame(x = runif(n1),
                      y = runif(n1))
coords2 <- data.frame(x = runif(n2),
                      y = runif(n2))
samples1 <- rep(1,n1)

eval_metrics1 <- c("CD")
eval_metrics2 <- c("CD", "DC")
eval_metrics3 <- c("CD", "DC", "ncol")

out1 <- wSIR:::wSIROptimisation(exprs_train = x1, 
                                coords_train = coords1, 
                                exprs_test = x2, 
                                coords_test = coords2, 
                                samples_train = samples1, 
                                slices = 5, 
                                alpha = 4, 
                                eval_metrics = eval_metrics1)
out2 <- wSIR:::wSIROptimisation(exprs_train = x1, 
                                coords_train = coords1, 
                                exprs_test = x2, 
                                coords_test = coords2, 
                                samples_train = samples1, 
                                slices = 5, 
                                alpha = 4, 
                                eval_metrics = eval_metrics2)
out3 <- wSIR:::wSIROptimisation(exprs_train = x1, 
                                coords_train = coords1, 
                                exprs_test = x2, 
                                coords_test = coords2, 
                                samples_train = samples1, 
                                slices = 5, 
                                alpha = 4, 
                                eval_metrics = eval_metrics3)
# check all right lengths
expect_equal(length(out1), 1)
expect_equal(length(out2), 2)
expect_equal(length(out3), 3)

