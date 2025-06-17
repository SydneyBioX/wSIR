context("Test createWeightMatrix")

nslices = 3
set.seed(12345)
coords = data.frame("x" = runif(n = 100, min = 0, max = 1),
                    "y" = runif(n = 100, min = 0, max = 1))
labels = spatialAllocator(coords = coords, slices = nslices)

wmat0 = createWeightMatrix(coords = coords, labels = labels, alpha = 0)
diag(wmat0) = 0

zero_matrix = matrix(rep(0, nslices^4), nrow = nslices^2)
expect_equal(wmat0, zero_matrix)