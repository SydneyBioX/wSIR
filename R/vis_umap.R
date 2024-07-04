#' vis_umap
#'
#' @description
#' A function to generate UMAP plots on the low-dimensional embedding of the gene expression data. The points are coloured by
#' their value for the genes with highest (in absolute value) loading in WSIR1.
#'
#' @param exprs matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param WSIR wsir object that is output of wSIR function. If you wish to generate UMAP plots based on other DR methods, ensure
#' that the slot named "scores" in WSIR parameter contains the low-dimensional representation of exprs.
#' @param highest_genes output from top_genes function.
#' @param n_genes integer for the number of genes you would like to show. Default is the number of genes you selected in the
#' top_genes function.
#'
#' @return Grid of umap plots with n_genes number of plots. Each shows the cells in a UMAP generated on the low-dimensional gene
#' expression data, coloured by their value for each of the genes found by top_genes.
#'
#' @importFrom umap umap
#'
#' @examples
#' wsir_obj = wSIR(exprs = sample1_exprs, coords = sample1_coords) # create wsir object
#' top_genes_obj = top_genes(WSIR = wsir_obj, highest = 4) # create top genes object
#' umap_plots = vis_umap(exprs = sample1_exprs, WSIR = wsir_obj, highest_genes = top_genes_obj, n_genes = 4)
#' umap_plots
#'
#' @export

vis_umap <- function(exprs,
                     WSIR,
                     highest_genes,
                     n_genes = length(names(highest_genes$genes))) {
  n_genes <- min(n_genes, length(names(highest_genes$genes))) # make sure it is a valid number of genes
  gene_names <- names(highest_genes$genes)[1:n_genes]

  umap_obj <- umap::umap(WSIR$scores) # create umap object
  gene_inds <- match(gene_names, colnames(exprs)) # identify the relevant columns in the gene expression matrix

  # create and fill umap_df
  umap_df <- matrix(NA, nrow = n_genes*nrow(exprs), ncol = 4) %>% as.data.frame()
  colnames(umap_df) <- c("UMAP1", "UMAP2", "gene", "expression")

  umap_df$UMAP1 <- rep(umap_obj$layout[,1], n_genes)
  umap_df$UMAP2 <- rep(umap_obj$layout[,2], n_genes)
  umap_df$gene <- vec_rep_each(gene_names, nrow(exprs))
  umap_df$expression <- exprs[, gene_inds] %>% as.vector()

  # create plot
  plot = ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, colour = expression)) +
    geom_point() +
    facet_wrap(~gene)
  return(plot)
}
