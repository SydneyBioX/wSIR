#context("Test generateUmapFromWSIR")

library(wSIR)
n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
colnames(x) <- paste0("g", c(1:p))
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

wsir_out <- wSIR(X = x,
                 coords = coords)
umap_coords <- generateUmapFromWSIR(wsir = wsir_out)

expect_equal(nrow(umap_coords), n)
expect_equal(ncol(umap_coords), 2)