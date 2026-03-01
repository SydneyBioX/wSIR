# spatialAllocator

This function allocates each cell to a tile based on a specified number
of tiles.

## Usage

``` r
spatialAllocator(coords, slices = 3)
```

## Arguments

- coords:

  dataframe contains the spatial position of each cell. Column names
  c("x", "y'). Must include row names as integer from 1 to nrow(coords).
  This is automatically included in the wSIR function prior to this
  function.

- slices:

  integer the number of slices along each of the two spatial axes. The
  number of tiles will be slices^2 since there are two spatial axes. E.g
  set slices = 4 to use 4^2 = 16 tiles. Default value is 3.

## Value

output by itself is a matrix containing slice belonging for each axis in
long format. When used in lapply(coords_split, spatialAllocator, slices)
as in wSIR the output is a dataframe with each cell's tile allocation in
the "coordinate" column.
