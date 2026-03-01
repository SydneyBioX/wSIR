# sirPCA

This function performs eigendecomposition on `(X^H)^t %*% W %*% (X^H)`,
where `X^H` is matrix of scaled slice means.

## Usage

``` r
sirPCA(
  sliced_data,
  maxDirections = 50,
  W = diag(nrow(sliced_data)),
  varThreshold = 0.99
)
```

## Arguments

- sliced_data:

  matrix of scaled slice means

- maxDirections:

  integer (default 50), number of directions we want in our final
  low-dimensional Z.

- W:

  matrix of slice weights. Output of createWeightMatrix function.

- varThreshold:

  numeric proportion of eigenvalues of variance in
  `t(X_H) %*% W %*% X_H` to retain. Default is 0.99. Select higher
  threshold to include more dimensions, lower threshold to include less
  dimensions.

## Value

list containing: 1) "eigenvectors" matrix of eigenvectors of
`(X^H)^t %*% W %*% (X^H)` 2) "d" integer, selected number of dimensions
3) "evalues" vector of eigenvalues of `(X^H)^t %*% W %*% (X^H)`
