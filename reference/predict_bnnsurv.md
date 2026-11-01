# Predict Survival with a Bagged Nearest-Neighbour Model

Generates survival probabilities from an `mlsurv_model` fitted by
[`fit_bnnsurv`](https://elbadisyimad.com/survalis/reference/fit_bnnsurv.md).

## Usage

``` r
predict_bnnsurv(object, newdata, times = NULL)
```

## Arguments

- object:

  A fitted `mlsurv_model` returned by
  [`fit_bnnsurv()`](https://elbadisyimad.com/survalis/reference/fit_bnnsurv.md).

- newdata:

  A data frame of new observations for prediction.

- times:

  Optional numeric vector of evaluation time points. If `NULL`
  (default), survival is returned at every distinct training time.

## Value

A base `data.frame` with one row per observation and columns named
`"t=<time>"` containing survival probabilities in `[0, 1]`,
non-increasing over time.

## Details

The ensemble curve is a step function that jumps only at training event
times, so the value at time \\t\\ is the curve at the last training time
not after \\t\\ (and `1` before the first one). Factor levels in
`newdata` are matched to the training levels; unseen levels and missing
predictor values are errors.

## See also

[`fit_bnnsurv()`](https://elbadisyimad.com/survalis/reference/fit_bnnsurv.md),
[`tune_bnnsurv()`](https://elbadisyimad.com/survalis/reference/tune_bnnsurv.md)

## Examples

``` r
mod <- fit_bnnsurv(Surv(time, status) ~ age + karno + diagtime + prior, data = veteran)
pred <- predict_bnnsurv(mod, newdata = veteran[1:3, ], times = c(50, 100, 200))
pred
#>   t=50 t=100 t=200
#> 1 0.80  0.34  0.24
#> 2 0.96  0.52  0.22
#> 3 0.82  0.40  0.40
```
