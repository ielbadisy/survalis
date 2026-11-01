# Predict Survival Probabilities from a Conditional Inference Survival Forest

Generates predicted survival probabilities at specified time points from
a fitted `cforest` survival model.

## Usage

``` r
predict_cforest(object, newdata, times, ...)
```

## Arguments

- object:

  An `mlsurv_model` object returned by
  [`fit_cforest()`](https://elbadisyimad.com/survalis/reference/fit_cforest.md).

- newdata:

  A data frame containing the predictor variables for prediction.

- times:

  A numeric vector of time points at which to estimate survival
  probabilities.

- ...:

  Not used, included for compatibility.

## Value

A data frame of survival probabilities with one row per observation in
`newdata` and one column per requested time point (named `"t=<time>"`).

## Details

Survival curves are extracted from the fitted `cforest` model and
linearly interpolated to the requested time points.

## Examples

``` r
mod <- fit_cforest(Surv(time, status) ~ age + celltype + karno, data = veteran)
predict_cforest(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#>       t=100     t=200     t=300
#> 1 0.5359829 0.2793382 0.1053863
#> 2 0.8196588 0.4229460 0.2956962
#> 3 0.5441606 0.4260393 0.1918215
#> 4 0.5948025 0.3320945 0.1473169
#> 5 0.8071864 0.3759406 0.2538033
```
