context("Test plotUmapFromWSIR")

library(wSIR)
n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
colnames(x) <- paste0("g", c(1:p))
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

wsir_out <- wSIR(X = x,
                 coords = coords)
umap_coords <- generateUmapFromWSIR(wsir = wsir_out)
top_genes_obj <- findTopGenes(wsir = wsir_out, highest = 4)

umap_plot <- plotUmapFromWSIR(umap_coords = umap_coords,
                              X = x,
                              highest_genes = top_genes_obj,
                              n_genes = 4)

# expect ggplot object. Also this tests if there is a plot itself 
expect_true("ggplot" %in% class(umap_plot)) 