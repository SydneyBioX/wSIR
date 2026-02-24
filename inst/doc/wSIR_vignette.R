## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = FALSE
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(wSIR) # package itself
library(magrittr) # for %>% 
library(ggplot2) # for ggplot
library(doBy) # for which.maxn
library(vctrs) # for vec_rep_each
library(umap) # for umap
library(class) # for example wSIR application

## -----------------------------------------------------------------------------
data(MouseData)

## -----------------------------------------------------------------------------
# a <- Sys.time()
# optim_obj <- exploreWSIRParams(X = as.matrix(sample1_exprs), 
#                                coords = sample1_coords, 
#                                # optim_alpha = c(0,.5, 1,2,4,8), 
#                                # optim_slices = c(5,10,15),
#                                optim_alpha = c(0, 4),
#                                optim_slices = c(10, 15),
#                                nrep = 5,
#                                metric = "DC")
# Sys.time()-a
# optim_obj$plot

## -----------------------------------------------------------------------------
wsir_obj <- wSIR(X = sample1_exprs, 
                 coords = sample1_coords, 
                 slices = 10,#optim_obj$best_slices, 
                 alpha = 4,#optim_obj$best_alpha,
                 optim_params = FALSE)

names(wsir_obj)

## -----------------------------------------------------------------------------
top_genes_obj <- findTopGenes(WSIR = wsir_obj, highest = 8) # create top genes object
top_genes_plot <- top_genes_obj$plot # select plot
top_genes_plot # print plot

top_genes_obj <- findTopGenes(WSIR = wsir_obj, highest = 8, dirs = 2:4)
top_genes_plot <- top_genes_obj$plot
top_genes_plot

## -----------------------------------------------------------------------------
vis_obj <- visualiseWSIRDirections(coords = sample1_coords, WSIR = wsir_obj, dirs = 8) # create visualisations
vis_obj

## -----------------------------------------------------------------------------
umap_coords <- generateUmapFromWSIR(WSIR = wsir_obj)
umap_plots <- plotUmapFromWSIR(X = sample1_exprs,
                               umap_coords = umap_coords,
                               highest_genes = top_genes_obj,
                               n_genes = 6)
umap_plots

## -----------------------------------------------------------------------------
sample2_low_dim_exprs <- projectWSIR(wsir = wsir_obj, newdata = sample2_exprs)

## -----------------------------------------------------------------------------
dim(sample2_low_dim_exprs)

## -----------------------------------------------------------------------------
head(sample2_low_dim_exprs)

## -----------------------------------------------------------------------------
wsir_obj_samples12 <- wSIR(X = rbind(sample1_exprs, sample2_exprs),
                           coords = rbind(sample1_coords, sample2_coords),
                           samples = c(rep(1, nrow(sample1_coords)), rep(2, nrow(sample2_coords))),
                           slices = 10,#optim_obj$best_slices, 
                           alpha = 4,#optim_obj$best_alpha,
                           optim_params = FALSE)

sample3_low_dim_exprs <- projectWSIR(wsir = wsir_obj_samples12, newdata = sample3_exprs)
dim(sample3_low_dim_exprs)

## -----------------------------------------------------------------------------
samples12_cell_types <- append(sample1_cell_types, sample2_cell_types)

knn_classification_object <- knn(train = wsir_obj_samples12$scores, 
                                 test = sample3_low_dim_exprs,
                                 cl = samples12_cell_types,
                                 k = 10)

tail(knn_classification_object)

## -----------------------------------------------------------------------------
sessionInfo()

