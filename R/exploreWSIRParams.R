# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param X no arguments
#' @param coords to fill
#' @param alphaVals to fill
#' @param numSlices to fill
#' @param nRepeats to fill
#' @param ... to fill
#'
#' @return prints hello world
#'
#' @examples
#' #hello()
#'
#' @export
exploreWSIRParams = function(X,
                             coords,
                             alphaVals = seq(from = -10, to = 10, by = 1),
                             numSlices = 3:10,
                             nRepeats = 1,
                             ...) {
  # this function performs repeated wSIR for different values of
  # slices and alpha and will output the results as well as diagnostic
  # plot to help the user decide which parameters should be chosen.
  # overall those with the highest distCor should be selected


  # ... passes onto wSIR function

  # example usage:
  # expl = exploreWSIRParams(X,
  #                          coords,
  #                          alphaVals = seq(from = -5, to = 5, by = 1),
  #                          numSlices = 3:10,
  #                          nRepeats = 5,
  #                          maxDirections = 20,
  #                          varThreshold = 0.9)

  # require(Rfast)
  # require(patchwork)

  # res_df_list = list()
  res_array = array(data = NA, dim = c(length(alphaVals), length(numSlices), nRepeats),
                    dimnames = list(alphaVals, numSlices, seq_len(nRepeats)))

  for (rep in seq_len(nRepeats)) {

    print(rep)

    # step 1 split the data X into a training and validation set

    keep = sample(c(TRUE, FALSE), nrow(coords), replace = TRUE)

    X_train = X[keep,]
    coords_train <- coords[keep,]

    X_val <- X[!keep,]
    coords_val <- coords[!keep,]

    for (i in seq_len(length(alphaVals))) {
      # print(i)
      for (j in seq_len(length(numSlices))) {
        # print(j)
        out = tryCatch(wSIR(X = X_train,
                            coords = coords_train,
                            slices = numSlices[j],
                            weighted = TRUE,
                            alpha = exp(alphaVals[i]),
                            ...), error = function(e) return(NULL))
        if (is.null(out)) next
        proj = project_SIR(out, as.matrix(X_val))
        res_array[i, j, rep] <- Rfast::dcor(x = coords_val, y = proj)$dcor
      }
    }

  }

  res_array_mean = apply(res_array, c(1,2), mean)

  res_df = as.data.frame.table(res_array_mean)
  colnames(res_df) <- c("logAlpha", "numSlices", "distCor")
  res_df$logAlpha <- as.numeric(as.character(res_df$logAlpha))
  res_df$numSlices <- as.numeric(as.character(res_df$numSlices))

  g1 = ggplot(res_df, aes(x = logAlpha, y = numSlices)) +
    geom_tile(aes(fill = distCor)) +
    scale_fill_gradient2(low = "grey", mid = "red", high = "black",
                         na.value = "white", midpoint = mean(res_df$distCor)) +
    theme_classic()

  g2 = ggplot(res_df, aes(x = logAlpha, y = distCor)) +
    geom_line(aes(group = numSlices, colour = numSlices), linetype = "solid") +
    scale_colour_gradient(low = "grey", high = "black", na.value = "white") +
    theme_classic()

  gBoth = g1 + g2

  return(list(res_array = res_array, res_df = res_df, plot = gBoth))

}
