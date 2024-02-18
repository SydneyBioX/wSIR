# Hello world function

#' Day hello
#'
#' @description
#' This function says hello
#'
#' @param none no arguments
#'
#' @return prints hello world
#'
#' @examples
#' hello()
#'
#' @export
slicer <- function(X, Y, slices = 8, categorical = FALSE) { # this is a slicer for univariate Y: if Y is multivariate, need to make it univariate in a previous step

  ## slicer

  ### inputs:
  # X: dataframe or matrix
  # Y: dataframe with one column
  # slices: integer with number of slices. (Only) necessary if response is continuous.
  # categorical: binary, states if the response is categorical or not. Default is FALSE ( = continuous response)

  ### output: dataframe of size slices * p, containing the means for each slice for each column.


  # Ensure X is a data.frane
  X <- X %>% as.data.frame()
  if (!categorical) { # if Y continuous, create new 'slice' column in X in the appropriate way
    data <- sorter(X, Y)
    nX <- nrow(X)
    data$slice <- allocator(nrows = nX, slices = slices)
  } else { # if Y categorical: add new 'slice' column that is Y values (and sort by those values)
    data <- X %>% as.data.frame()
    data$slice <- Y[,1]
    data <- data %>% arrange(slice)
    slices <- length(unique(data$slice))
  }
  slice_names <- unique(data$slice) # need this change so that it works when we have slice names as tile coords

  long_values_sliced_dataframe <- c()

  for (s in 1:slices) {
    for (i in 1:nrow(data)) {
      if (data$slice[i] == slice_names[s]) {
        first_row <- i
        break
      }
    }
    for (i in 1:nrow(data)) {
      if (data$slice[i] == slice_names[s]) {
        last_row <- i
      }
    }
    dataset_for_avging <- data[first_row:last_row,1:ncol(X)]
    avgs <- means(dataset_for_avging)
    long_values_sliced_dataframe <- long_values_sliced_dataframe %>% append(avgs)
  }
  sliced_dataset <- matrix(long_values_sliced_dataframe, ncol = ncol(X), byrow = TRUE) %>% as.data.frame()
  #rownames(sliced_dataset) <- c(1:slices)
  return(sliced_dataset)
}
