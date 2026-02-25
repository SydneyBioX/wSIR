context("Test spatialAllocator")

library(wSIR)
# create extreme example of points being not together, i.e in a ring shape

slices = 6 # number of slices
n <- 200
minr <- 2 # inner radius of ring
maxr <- 3 # outer radius of ring
thetas <- runif(n = n, min = 0, max = 2*pi)
radii <- runif(n = n, min = minr, max = maxr)

xcoords <- radii * cos(thetas)
ycoords <- radii * sin(thetas)

coords <- data.frame("x" = xcoords,
                     "y" = ycoords)
rownames(coords) <- c(1:n)

allocation <- wSIR:::spatialAllocator(coords = coords, slices = slices)
n_tiles <- allocation$coordinate %>% unique() %>% length() # find number of unique tiles
# finish with an expect less tiles than 36 statement 