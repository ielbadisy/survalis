# Predict Survival Probabilities from a rpsurv Model

Generates predicted survival probabilities at the requested `times` for
a model fitted with
[`fit_rpsurv`](https://elbadisyimad.com/survalis/reference/fit_rpsurv.md).

## Usage

``` r
predict_rpsurv(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` from
  [`fit_rpsurv`](https://elbadisyimad.com/survalis/reference/fit_rpsurv.md).

- newdata:

  A `data.frame` of predictors (raw covariate columns, not pre-expanded
  dummies).

- times:

  Numeric vector of evaluation time points.

## Value

A `data.frame` of survival probabilities with one row per observation in
`newdata` and one column per time point (columns named `"t=<time>"`).

## See also

[`fit_rpsurv`](https://elbadisyimad.com/survalis/reference/fit_rpsurv.md)

## Examples

``` r
mod <- fit_rpsurv(Surv(time, status) ~ age + karno, veteran, df = 3)
predict_rpsurv(mod, veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200      t=300
#> 1 0.4012090 0.1792489 0.08348534
#> 2 0.5193893 0.2914039 0.16844492
#> 3 0.3910155 0.1707731 0.07784344
#> 4 0.3992378 0.1775949 0.08237479
#> 5 0.5196942 0.2917260 0.16871396
```
