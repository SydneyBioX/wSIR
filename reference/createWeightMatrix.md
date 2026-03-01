# createWeightMatrix

A function to create the weight matrix given the location of the cells,
tile allocations and desired spatial weighting strength. Weight matrix
entries represent level of spatial correlation between all pairs of
tiles.

## Usage

``` r
createWeightMatrix(coords, labels, alpha = 4)
```

## Arguments

- coords:

  dataframe of dimension n \* 2. Column names c("x", "y"). Spatial
  position of each cell.

- labels:

  dataframe of dimension n \* 1, column name c("coordinate"). Tile
  allocation of each cell. This is automatically created in the wSIR
  function.

- alpha:

  numeric value giving strength of spatial weighting matrix. alpha = 0
  returns identity matrix and equals SIR. Large alpha values tend all
  entries towards 1. Default is 4.

## Value

matrix containing the weight value for all pairs of tiles. Each value is
between 0 and 1, with 1 always on the diagonal.
