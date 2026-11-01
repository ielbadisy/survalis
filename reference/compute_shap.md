# Compute local SHAP-like contributions for survival predictions

Estimates per-feature contributions (à la SHAP) to the predicted
survival probability for a **single** observation at one or more time
points. For each time, the method samples random feature orders,
marginalizes future features with values from a baseline dataset, and
accumulates the marginal effects of adding each feature back.

## Usage

``` r
compute_shap(
  model,
  newdata,
  baseline_data,
  times,
  sample.size = 100,
  features = NULL,
  aggregate = FALSE,
  method = c("meanabs", "integral")
)
```

## Arguments

- model:

  A fitted survival model produced by `fit_*()` that exposes
  `model$learner` (e.g., "coxph") and `model$data` (the training frame).

- newdata:

  A data frame with **exactly one row** (the instance to explain).

- baseline_data:

  A data frame to sample background values from (typically the training
  data used to fit `model`).

- times:

  Numeric vector of evaluation times (same scale as the outcome).

- sample.size:

  Integer, number of random feature orderings to sample per time
  (default `100`).

- features:

  Optional character vector of covariate names to explain. If `NULL`
  (default), inferred as
  `setdiff(colnames(newdata), c(model$time, model$status))`, consistent
  with
  [`compute_interactions`](https://elbadisyimad.com/survalis/reference/compute_interactions.md).

- aggregate:

  Logical; if `TRUE`, aggregate contributions across `times` into a
  single value per feature (see `method`).

- method:

  Character; aggregation method if `aggregate=TRUE`: `"meanabs"` (mean
  absolute contribution) or `"integral"` (trapezoidal integral over
  time). Default `"meanabs"`.

## Value

If `aggregate = FALSE`: a data frame with columns `feature`, `phi`, and
`time` (one row per feature per time). If `aggregate = TRUE`: a data
frame with columns `feature` and `phi` (one row per feature), with
attribute `"shap_method"` set to the aggregation used.

## Details

The prediction function is inferred from `model$learner` as
`predict_<learner>` and called with signature
`predict_fun(model, newdata, times = times)`. Factor levels in `newdata`
are harmonized to those in `model$data`.

## See also

[`plot_shap()`](https://elbadisyimad.com/survalis/reference/plot_shap.md),
the various `predict_*` methods (e.g.
[`predict_coxph()`](https://elbadisyimad.com/survalis/reference/predict_coxph.md))

## Examples

``` r
# \donttest{
mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)

shap_td <- compute_shap(
  model         = mod,
  newdata       = veteran[100, , drop = FALSE],
  baseline_data = veteran,
  times         = c(100, 200),
  sample.size   = 5,
  aggregate     = FALSE
)
head(shap_td)
#>        feature         phi time
#> 100.1      trt  0.00000000  100
#> 100.2 celltype -0.02899675  100
#> 100.3    karno  0.00000000  100
#> 100.4 diagtime  0.00000000  100
#> 100.5      age  0.00000000  100
#> 100.6    prior  0.00000000  100

shap_meanabs <- compute_shap(
  model         = mod,
  newdata       = veteran[100, , drop = FALSE],
  baseline_data = veteran,
  times         = c(100, 200),
  sample.size   = 5,
  aggregate     = TRUE,
  method        = "meanabs"
)
head(shap_meanabs)
#>           feature         phi
#> trt           trt 0.000000000
#> celltype celltype 0.023253617
#> karno       karno 0.035206027
#> diagtime diagtime 0.000000000
#> age           age 0.007297201
#> prior       prior 0.000000000
# }
```
