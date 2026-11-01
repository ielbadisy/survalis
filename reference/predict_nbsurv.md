# Predict Survival Probabilities from a nbsurv Model

Generates predicted survival probabilities at the requested `times` for
a model fitted with
[`fit_nbsurv`](https://elbadisyimad.com/survalis/reference/fit_nbsurv.md).

## Usage

``` r
predict_nbsurv(object, newdata, times)
```

## Arguments

- object:

  An `"mlsurv_model"` from
  [`fit_nbsurv`](https://elbadisyimad.com/survalis/reference/fit_nbsurv.md).

- newdata:

  A `data.frame` of predictors.

- times:

  Numeric vector of evaluation time points.

## Value

A `data.frame` of survival probabilities with one row per observation in
`newdata` and one column per time point (columns named `"t=<time>"`).

## See also

[`fit_nbsurv`](https://elbadisyimad.com/survalis/reference/fit_nbsurv.md)

## Examples

``` r
mod <- fit_nbsurv(Surv(time, status) ~ age + karno, veteran)
predict_nbsurv(mod, veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200      t=300
#> 1 0.4190910 0.1591754 0.08977703
#> 2 0.6179400 0.2755180 0.17334416
#> 3 0.4621632 0.3020846 0.08723375
#> 4 0.6203305 0.3152120 0.19393861
#> 5 0.6213628 0.2716229 0.17120548
```
