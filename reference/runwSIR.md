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

  arguments passing to `calculateWSIR`

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
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
library(SpatialExperiment)

sce <- SingleCellExperiment(assays = list(logcounts = t(sample1_exprs)),
reducedDims = list(spatial = sample1_coords))

sce <- runwSIR(x = sce, dim_red = "spatial")

spe <- SpatialExperiment(assays = list(logcounts = t(sample1_exprs)),
spatialCoords = as.matrix(sample1_coords))

spe <- runwSIR(x = spe, spatialCoords = TRUE)
```
