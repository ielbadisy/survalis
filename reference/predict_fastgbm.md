# Predict Survival Probabilities from a fastgbm Model

Generates predicted survival probabilities at the requested `times` for
a model fitted with
[`fit_fastgbm`](https://elbadisyimad.com/survalis/reference/fit_fastgbm.md)
(piecewise-exponential objective).

## Usage

``` r
predict_fastgbm(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` from
  [`fit_fastgbm`](https://elbadisyimad.com/survalis/reference/fit_fastgbm.md).

- newdata:

  A `data.frame` of predictors.

- times:

  Numeric vector of evaluation time points.

## Value

A `data.frame` of survival probabilities with one row per observation in
`newdata` and one column per time point (columns named `"t=<time>"`).

## See also

[`fit_fastgbm`](https://elbadisyimad.com/survalis/reference/fit_fastgbm.md)

## Examples

``` r
mod <- fit_fastgbm(Surv(time, status) ~ age + karno, veteran, ntrees = 50)
predict_fastgbm(mod, veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.6239460 0.4099494 0.2789265
#> 2 0.7262066 0.5117303 0.3983366
#> 3 0.6545936 0.4101674 0.2798156
#> 4 0.6492212 0.4556793 0.3320398
#> 5 0.7206470 0.4954610 0.3789357
```
