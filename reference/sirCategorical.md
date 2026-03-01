# sirCategorical

This function performs WSIR based on provided gene expression data, tile
allocations, and weight matrix.

## Usage

``` r
sirCategorical(X, Y, W = NULL, ...)
```

## Arguments

- X:

  matrix of normalised gene expression data for n genes across p cells.

- Y:

  dataframe with 1 column named "coordinate" that is the tile allocation
  for each cell. There should be up to slices^2 unique tile IDs in this
  column.

- W:

  Weight matrix created by createWeightMatrix. Entry (i,j) represents
  the spatial correlation level between tiles i and j. The diagonal
  values should be all 1. If not provided, SIR implementation will be
  used.

- ...:

  arguments passed to sirPCA

## Value

list of outputs with 5 named slots. They are the same as the output of
the wSIR function: this is the final step in the wSIR function.
