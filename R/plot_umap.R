#' plot_umap
#'
#' @description
#' A function to plot a UMAP generated on the low-dimensional embedding of the gene expression data. The points are coloured by
#' their value for the genes with highest (in absolute value) loading in a selected WSIR direction, by default WSIR1.
#'
#' @param exprs matrix containing normalised gene expression data including n cells and p genes, dimension n * p.
#' @param umap_coords UMAP coordinates for each cell that is output of generate_umap function. The UMAP coordinates can be
#' based on any dimension reduction method, e.g they could be the UMAP coordinates computed on the WSIR dimension reduction
#' of the gene expression data, or on the PCs (principal components), or on any other low-dimensional matrix. Must be
#' a matrix of dimension nrow(exprs) * 2.
#' @param highest_genes output from top_genes function. Ensure that you only used a single WSIR direction in the top_genes
#' function, otherwise some genes are likely to be displayed multiple times in the resulting plots from this function.
#' Default is NULL so an error message can easily be thrown if genes and highest_genes are both not provided.
#' @param genes vector with gene names (must all be in colnames(exprs)) you wish to show in the UMAP plot. The cells
#' in those plots will be coloured by their expression values for the genes you provide here. Must provide either genes
#' or highest_genes parameter (not both): provide genes if you want to visualise a few specific genes, provide highest_genes
#' if you want to visualise the genes that are found to be the most important to the WSIR directions. Default is NULL
#' so an error message can easily be thrown if genes and highest_genes are both not provided.
#' @param n_genes integer for the number of genes you would like to show. Default is the number of unique genes in the
#' highest_genes parameter or the number of genes in the (vector) parameter genes. Use this parameter if you want to
#' show only a few of the most important genes (e.g select the top 4 with n_genes = 4).
#' @param ... additional parameters for ggplot functions, e.g size (for size of points).
#'
#' @return Grid of umap plots with n_genes number of plots. Each shows the cells in a UMAP generated on the low-dimensional gene
#' expression data, coloured by their value for each of the genes found by top_genes.
#'
#' @importFrom vctrs vec_rep_each
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom magrittr %>%
#'
#' @examples
#' data(MouseData)
#' wsir_obj = wSIR(X = sample1_exprs,
#'   coords = sample1_coords,
#'   optim_params = FALSE,
#'   alpha = 4,
#'   slices = 6) # create wsir object
#' umap_coords = generate_umap(WSIR = wsir_obj)
#' top_genes_obj = top_genes(WSIR = wsir_obj, highest = 4) # create top genes object
#' umap_plot = plot_umap(umap_coords = umap_coords,
#'   exprs = m1_exprs,
#'   highest_genes = top_genes_obj,
#'   n_genes = 4)
#' umap_plot
#'
#' @export

plot_umap <- function(exprs,
                     umap_coords,
                     highest_genes = NULL,
                     genes = NULL,
                     n_genes,
                     ...) {
  if (is.null(highest_genes) == is.null(genes)) { # error message if incorrect inputs
    return("Must provide one of highest_genes or genes, not neither nor both.")
  }

  if (is.null(genes)) {
    genes = unique(highest_genes$genes$gene) # if highest_genes is provided, extract the top genes from that parameter
  }

  n_genes <- min(n_genes, length(genes)) # make sure it is a valid number of genes
  gene_names <- genes[1:n_genes]

  gene_inds <- match(gene_names, colnames(exprs)) # identify the relevant columns in the gene expression matrix

  # create and fill umap_df
  umap_df <- matrix(NA, nrow = n_genes*nrow(exprs), ncol = 4) %>% as.data.frame()
  colnames(umap_df) <- c("UMAP1", "UMAP2", "gene", "expression")

  umap_df$UMAP1 <- rep(umap_coords[,1], n_genes)
  umap_df$UMAP2 <- rep(umap_coords[,2], n_genes)
  umap_df$gene <- vec_rep_each(gene_names, nrow(exprs))
  umap_df$expression <- as.matrix(exprs)[, gene_inds] %>% as.vector()

  # create plot
  plot = ggplot(data = umap_df, aes(x = UMAP1, y = UMAP2, colour = expression)) +
    geom_point() +
    facet_wrap(~gene) +
    theme_classic()
  return(plot)
}
