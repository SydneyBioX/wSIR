# runwSIR

Perform wSIR on cells, based on the expression data and a reducedDim in
a SingleCellExperiment or SpatialExperiment object

## Usage

``` r
runwSIR(x, name = "wSIR", scores_only = FALSE, ...)
```

## Arguments

- x:

  A numeric matrix of normalised gene expression data where rows are
  features and columns are cells. Alternatively, a SingleCellExperiment
  or SpatialExperiment containing such a matrix

- name:

  string to specify the name to store the result in the reducedDims of
  the output. Default is "wSIR"

- scores_only:

  logical whether only the wSIR scores should be calculated. If FALSE
  additional information about the wSIR model will be stored in the
  attributes of the object. Default FALSE.

- ...:

  arguments passing to `calculatewSIR`

## Value

If `x` is matrix-like, a list containing wSIR scores, loadings, etc. If
`x` is a SingleCellExperiment or SpatialExperiment, the same object is
returned with an additional slot in `reducedDims(..., name)`
corresponding to the wSIR scores matrix. If `scores_only = FALSE`, then
the attributes of the wSIR scores contain the following elements:

- directions

- estd

- W

- evalues

## Examples

``` r
data(MouseData)
library(SingleCellExperiment)
library(SpatialExperiment)

sce <- SingleCellExperiment(assays = list(logcounts = t(sample1_exprs)),
reducedDims = list(spatial = sample1_coords))

sce <- runwSIR(x = sce, dimred = "spatial")

spe <- SpatialExperiment(assays = list(logcounts = t(sample1_exprs)),
spatialCoords = as.matrix(sample1_coords))

spe <- runwSIR(x = spe, spatialCoords = TRUE)
```
