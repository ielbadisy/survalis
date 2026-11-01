# Predict Survival with XGBoost

Generates predictions from an XGBoost survival `mlsurv_model`. If
`times` is `NULL`, returns the raw linear predictor. If `times` is
provided, returns survival probabilities computed via the AFT mapping
using `object$dist` and `object$scale`.

## Usage

``` r
predict_xgboost(object, newdata, times = NULL)
```

## Arguments

- object:

  A fitted `mlsurv_model` returned by
  [`fit_xgboost()`](https://elbadisyimad.com/survalis/reference/fit_xgboost.md).

- newdata:

  A `data.frame` of new observations for prediction.

- times:

  Optional numeric vector of evaluation time points (same scale as
  training time). If `NULL`, returns the linear predictor.

## Value

- If `times` is `NULL`: a numeric vector of linear predictors (one per
  row of `newdata`).

- If `times` is provided: a base `data.frame` of survival probabilities
  with columns named `"t=<time>"` and row names `"ID_<i>"`.

## Details

AFT mapping with \\q = (\log t - \eta)/\sigma\\:

- Normal: \\S(t) = 1 - \Phi(q)\\

- Extreme (Gumbel): \\S(t) = \exp\\-\exp(q)\\\\

- Logistic: \\S(t) = 1/(1 + \exp(q))\\

Note: If the model was trained with `objective = "survival:cox"`,
survival probabilities are still computed using the supplied AFT
distribution/scale stored in the object; interpret with caution.

## See also

[`fit_xgboost()`](https://elbadisyimad.com/survalis/reference/fit_xgboost.md),
[`tune_xgboost()`](https://elbadisyimad.com/survalis/reference/tune_xgboost.md)

## Examples

``` r
# \donttest{
  mod_xgb <- fit_xgboost(Surv(time, status) ~ age + karno + celltype,
                         data = veteran, nrounds = 20)
  predict_xgboost(mod_xgb, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100      t=200      t=300
#> 1 0.3543772 0.12558322 0.04450383
#> 2 0.4031287 0.16251278 0.06551357
#> 3 0.3055987 0.09339059 0.02854005
#> 4 0.4008246 0.16066038 0.06439663
#> 5 0.4031287 0.16251278 0.06551357
# }
```
