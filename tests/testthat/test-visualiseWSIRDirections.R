context("Test visualiseWSIRDirections")

library(wSIR)
n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
colnames(x) <- paste0("g", c(1:p))
coords <- data.frame(x = runif(n),
                     y = runif(n))

wsir_out <- wSIR(X = x,
                 coords = coords)
vis_obj <- visualiseWSIRDirections(coords = coords,
                                   wsir = wsir_out, 
                                   dirs = 8) # create visualisations

# expect ggplot object. Also this tests if there is a plot itself 
expect_true("ggplot" %in% class(vis_obj)) 