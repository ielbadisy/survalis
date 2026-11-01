# Predict Survival Probabilities from a blackboost Model

Generates survival probabilities at requested times from a fitted mboost
Cox boosting model produced by
[`fit_blackboost()`](https://elbadisyimad.com/survalis/reference/fit_blackboost.md).

## Usage

``` r
predict_blackboost(object, newdata, times, ...)
```

## Arguments

- object:

  An `mlsurv_model` returned by
  [`fit_blackboost()`](https://elbadisyimad.com/survalis/reference/fit_blackboost.md).

- newdata:

  A `data.frame` of new observations for prediction.

- times:

  Numeric vector of time points at which to compute survival
  probabilities.

- ...:

  Additional arguments forwarded to
  [`mboost::survFit()`](https://rdrr.io/pkg/mboost/man/survFit.html).

## Value

A base `data.frame` with one row per observation in `newdata` and
columns named `"t=<time>"` (character), containing survival
probabilities in \[0, 1\].

## Details

Predictions use
[`mboost::survFit()`](https://rdrr.io/pkg/mboost/man/survFit.html) to
obtain a step function \\S(t)\\ per observation. If any requested
`times` exceed the model's maximum time, the last survival value is
carried forward (right-constant). Values are then matched to `times`
using stepwise (piecewise-constant) interpolation.

## Examples

``` r
mod <- fit_blackboost(Surv(time, status) ~ age + karno + celltype, data = veteran)
predict_blackboost(mod, newdata = veteran[1:5, ], times = c(5, 10, 40))
#>         t=5      t=10      t=40
#> 1 0.9750747 0.9287634 0.7361742
#> 2 0.9785472 0.9384807 0.7686258
#> 3 0.9750747 0.9287634 0.7361742
#> 4 0.9750747 0.9287634 0.7361742
#> 5 0.9785472 0.9384807 0.7686258
```
