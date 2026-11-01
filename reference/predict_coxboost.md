# Predict Survival with CoxBoost

Generates survival-probability predictions from a CoxBoost
`mlsurv_model`.

## Usage

``` r
predict_coxboost(object, newdata, times = NULL)
```

## Arguments

- object:

  A fitted `mlsurv_model` returned by
  [`fit_coxboost()`](https://elbadisyimad.com/survalis/reference/fit_coxboost.md).

- newdata:

  A `data.frame` of new observations for prediction.

- times:

  Optional numeric vector of evaluation time points. If `NULL`, returns
  the linear predictor.

## Value

- If `times` is `NULL`: a numeric vector of linear predictors (one per
  row of `newdata`).

- If `times` is provided: a base `data.frame` of survival probabilities
  with columns named `"t=<time>"`.

## See also

[`fit_coxboost()`](https://elbadisyimad.com/survalis/reference/fit_coxboost.md),
[`tune_coxboost()`](https://elbadisyimad.com/survalis/reference/tune_coxboost.md)

## Examples

``` r
# \donttest{
  mod_cb <- fit_coxboost(Surv(time, status) ~ age + karno + celltype,
                         data = veteran, stepno = 30)
  predict_coxboost(mod_cb, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.4422265 0.2030143 0.1087308
#> 2 0.5260552 0.2850004 0.1743259
#> 3 0.4422265 0.2030143 0.1087308
#> 4 0.4422265 0.2030143 0.1087308
#> 5 0.5260552 0.2850004 0.1743259
# }
```
