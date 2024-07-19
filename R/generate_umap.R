#' generate_umap
#'
#' @description
#' A function to generate UMAP coordinates from the low-dimensional embedding of the gene expression data. These
#' coordinates can later be plotted with the plot_umap function. Those two functions are separate so that you
#' can generate the UMAP points only once (which takes a long time), then modify the resulting plot as much as
#' desired.
#'
#' @param WSIR wsir object that is output of wSIR function. If you wish to generate UMAP plots based on other DR methods, ensure
#' that the slot named "scores" in WSIR parameter contains the low-dimensional representation of exprs.
#'
#' @return matrix of UMAP coordinates of dimension nrow(coords) * 2. Output of this function can be directly used
#' as the input to the plot_umap function.
#'
#' @importFrom umap umap
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
#'   exprs = sample1_exprs,
#'   highest_genes = top_genes_obj,
#'   n_genes = 4)
#' umap_plot
#'
#' @export

generate_umap <- function(WSIR) {
  umap_obj <- umap::umap(WSIR$scores) # create umap object
  return(umap_obj$layout)
}
