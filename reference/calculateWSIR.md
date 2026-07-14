# calculateWSIR

Perform wSIR on cells, based on the expression data and a reducedDim in
a SingleCellExperiment or SpatialExperiment object

## Usage

``` r
calculateWSIR(
  x,
  assay_type = "logcounts",
  dim_red = NULL,
  colData_columns = NULL,
  spatialCoords = FALSE,
  ...
)
```

## Arguments

- x:

  A numeric matrix of normalised gene expression data where rows are
  features and columns are cells. Alternatively, a SingleCellExperiment
  or SpatialExperiment containing such a matrix

- assay_type:

  if `x` is a SingleCellExperiment of SpatialExperiment then this is the
  assay for which wSIR will be calculated. Default "logcounts".

- dim_red:

  String or integer scalar specifying the dimensionality reduction slot
  for which to use for the slicing mechanism. Ignored if `coords` given.

- colData_columns:

  character vector specifying the subset of colData columns to be used
  for the wSIR slicing mechanism. Ignored if `coords` or `dim_red` given

- spatialCoords:

  logical indicating if spatialCoords should be used for the wSIR
  slicing mechanism. Ignored if `coords`, `dim_red`, or
  `colData_columns` given, or if `x` is not a SpatialExperiment object.

- ...:

  arguments passing to `wSIR`

## Value

A wSIR object
