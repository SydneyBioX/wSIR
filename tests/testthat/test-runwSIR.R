context("Test runwSIR")

library(wSIR)
library(SpatialExperiment)
library(SingleCellExperiment)

# Here test runwSIR function

## Test for a correctly named slot appearing in the output and that it is filled
## (will check that by checking it is correctly dimensioned)

n <- 100
p <- 10
x <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
coords <- matrix(runif(n = n*2), nrow = n, ncol = 2)

spe <- SpatialExperiment(
  assays = list(logcounts = t(x)),
  spatialCoords = coords
)

out <- runwSIR(x = spe,
               name = "wSIR_slot")

extract_lowdim <- reducedDim(out, "wSIR_slot")
expect_equal(nrow(extract_lowdim), n) # this checks that the reduced dim exists and is correct

## Now test for if the spatialCoords is missing
sce_missing <- SingleCellExperiment(assays = list(logcounts = t(x)))
expect_error(out <- runwSIR(x = sce_missing)) # get error if coords not present
