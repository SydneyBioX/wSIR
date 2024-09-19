#' visualiseWSIRDirections
#'
#' @description
#' A function to easily visualise the low-dimensional gene expression data. This function plots
#' each cell at its true spatial coordinates, coloured by their values for WSIR1 / WSIR2 / ... .
#' The plots give an intuition about what biological signals are contained in the WSIR directions.
#'
#' @param coords dataframe containing spatial positions of n cells in 2D space. Dimension n * 2. Column names must be c("x", "y").
#' @param WSIR wsir object as output of wSIR function. To analyse a different DR method, ensure the
#' slot named 'directions' contains the loadings as a matrix with the gene names as the rownames. Must
#' have used the same coords parameter as in coords parameter for this function.
#' @param dirs integer for how many of the low-dimensional directions you would like to visualise. Recommend no more
#' than 10 for ease of visualisation. Default is 6.
#' @param mincol String for the colour of low values of low-dimensional directions. Personal choice for user, default is "blue".
#' @param maxcol String for the colour of high values of low-dimensional directions. Personal choice for user, default is "red".
#'
#' @return Grid of plots with dirs number of plots. Each shows the cells at their spatial positions
#' coloured by their value for each of the first 'dirs' WSIR directions.
#'
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_point theme_classic facet_wrap ggtitle scale_color_gradient
#' @importFrom vctrs vec_rep_each
#' @importFrom rlang .data
#'
#' @examples
#' data(MouseData)
#' wsir_obj = wSIR(X = sample1_exprs,
#'   coords = sample1_coords,
#'   optim_params = FALSE,
#'   alpha = 4,
#'   slices = 6) # create wsir object
#' vis_obj = visualiseWSIRDirections(coords = sample1_coords,
#' WSIR = wsir_obj, dirs = 8) # create visualisations
#' vis_obj
#'
#' @export

visualiseWSIRDirections <- function(coords, WSIR, dirs = 6, mincol = "blue", maxcol = "red") {
  dirs <- min(dirs, ncol(WSIR$scores)) # make sure it is a valid value

  # initiliase empty long df
  vis_df_long <- matrix(NA, nrow = dirs*nrow(coords), ncol = 4) %>% as.data.frame()
  colnames(vis_df_long) <- c("x", "y", "value", "WSIR_direction")

  # fill columns of long df with relevant WSIR1/2/... values
  vis_df_long$x <- rep(coords$x, dirs)
  vis_df_long$y <- rep(coords$y, dirs)
  vis_df_long$value <- as.vector(WSIR$scores[,1:dirs])
  vis_df_long$WSIR_direction <- as.factor(vec_rep_each(c(1:dirs), nrow(coords)))

  # produce plot
  plot <- ggplot2::ggplot(aes(x = .data$x, y = .data$y, color = .data$value), data = vis_df_long) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::facet_wrap(~WSIR_direction, scales = "fixed") + # one panel per WSIR direction
    ggplot2::ggtitle("Cells at true positions coloured by WSIR values") +
    ggplot2::scale_color_gradient(low = mincol, high = maxcol)
  return(plot)
}
