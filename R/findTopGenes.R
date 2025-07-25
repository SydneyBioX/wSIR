#' findTopGenes
#'
#' @description
#' A function to find and visualise the genes with the highest
#' (in absolute value) loading in WSIR1.
#' These genes contribute the most to the first low-dimensional direction.
#'
#' @param WSIR wsir object as output of wSIR function. To analyse a different
#' DR method, ensure the
#' slot named 'directions' contains the loadings as a matrix with the gene
#' names as the rownames.
#' @param highest integer for how many of the top genes you would like to see.
#' Recommend no more
#' than 20 for ease of visualisation. Default is 10.
#' @param dirs integer or vector for which direction / directions you want to
#' show/find the top genes from.
#'
#' @return List containing two slots. First is named "plot" and shows a
#' barplot with the top genes
#' and their corresponding loadings. Second is named "genes" and is a vector
#' of the genes with highest
#' (in absolute value) loading in low-dimensional direction 1. Length is
#' parameter highest.
#'
#' @importFrom doBy which.maxn
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot aes geom_col theme_minimal labs geom_hline
#' @importFrom ggplot2 ggtitle facet_wrap
#' @importFrom stats reorder
#' @importFrom vctrs vec_rep_each
#' @importFrom rlang .data
#'
#' @examples
#' data(MouseData)
#'
#' wsir_obj = wSIR(X = sample1_exprs,
#'   coords = sample1_coords,
#'   optim_params = FALSE,
#'   alpha = 4,
#'   slices = 6) # create wsir object
#' top_genes_obj = findTopGenes(WSIR = wsir_obj, highest = 8)
#' top_genes_plot = top_genes_obj$plot # select plot
#' top_genes_plot # print plot
#'
#' @export

findTopGenes <- function(WSIR, highest = 10, dirs = 1) {

    wsir_dirs_df <- WSIR$directions %>%
        as.data.frame()
    wsir_dirs_df$gene <- rownames(WSIR$directions)

    res_df <- matrix(NA, nrow = length(dirs)*highest, ncol = 3) %>%
        as.data.frame()
    colnames(res_df) <- c("gene", "loading", "direction")
    res_df$direction <- vctrs::vec_rep_each(paste0("WSIR",dirs), highest)

    j <- 0
    for (i in dirs) {
        current_top_n <- doBy::which.maxn(abs(wsir_dirs_df[,i]), highest)
        current_genes <- wsir_dirs_df$gene[current_top_n]
        current_loadings <- wsir_dirs_df[current_top_n,i]
        res_df$gene[(j*highest+1):((j+1)*highest)] <- current_genes
        res_df$loading[(j*highest+1):((j+1)*highest)] <- current_loadings
        j <- j + 1
    }

    loadings_plot <- ggplot2::ggplot(
        aes(x = stats::reorder(.data$gene, -.data$loading, sum),
            y = .data$loading), data = res_df) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Gene", y = "Loading values from WSIR directions") +
        ggplot2::geom_hline(yintercept = 0) +
        ggplot2::ggtitle(paste0("Top ",
                                highest,
                                " genes with highest/lowest loading in wSIR ",
                                dirs)) +
        ggplot2::facet_wrap(~direction, nrow = 2, scales = "free")

    return(list(plot = loadings_plot,
                genes = res_df))
}
