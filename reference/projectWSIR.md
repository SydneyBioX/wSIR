# projectWSIR

function to project new gene expression data into low-dimensional space

## Usage

``` r
projectWSIR(wsir, newdata)
```

## Arguments

- wsir:

  wsir object that is usually the output of wSIR function. If you want
  to project new data into low-dim space following a different DR
  method, at param wsir use a list with matrix of loadings in slot 2
  (e.g PCA loadings) of dimension p \* d

- newdata:

  matrix of new gene expression data to project into low-dimensional
  space. Must have the same p columns as the columns in X argument used
  to generate wsir.

## Value

matrix of low-dimensional representation of newdata gene expression data

## Examples

``` r
data(MouseData)
wsir_obj <- wSIR(X = sample1_exprs,
    coords = sample1_coords,
    optim_params = FALSE,
    alpha = 4,
    slices = 6)
sample2_low_dim_exprs <- projectWSIR(wsir = wsir_obj, newdata = sample2_exprs)
```
