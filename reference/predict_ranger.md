# Predict Survival Probabilities from a ranger Model

Generates predicted survival probabilities for given time points from a
model fitted with
[`fit_ranger`](https://elbadisyimad.com/survalis/reference/fit_ranger.md).

## Usage

``` r
predict_ranger(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` object returned by
  [`fit_ranger`](https://elbadisyimad.com/survalis/reference/fit_ranger.md).

- newdata:

  A `data.frame` containing the same predictors as the training data.

- times:

  Numeric vector of time points at which to estimate survival
  probabilities.

## Value

A `data.frame` with one row per observation in `newdata` and one column
per time point (columns named `"t=<time>"`).

## Details

Predictions are obtained from
[`predict.ranger`](http://imbs-hl.github.io/ranger/reference/predict.ranger.md)
and survival curves are interpolated to match the requested `times`.

## See also

[`fit_ranger`](https://elbadisyimad.com/survalis/reference/fit_ranger.md),
[`tune_ranger`](https://elbadisyimad.com/survalis/reference/tune_ranger.md)

## Examples

``` r
mod <- fit_ranger(
  Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  num.trees = 25
)
predict_ranger(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.4319984 0.2594131 0.1846426
#> 2 0.7972877 0.4055384 0.3834524
#> 3 0.5123913 0.3480764 0.2339804
#> 4 0.7331676 0.3410541 0.2343888
#> 5 0.8619658 0.3637799 0.3495159
```
