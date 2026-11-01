# Predict Survival Probabilities from an RSF Model

Generates survival probability predictions for new data using a model
fitted by
[`fit_rsf()`](https://elbadisyimad.com/survalis/reference/fit_rsf.md).

## Usage

``` r
predict_rsf(object, newdata, times = NULL, ...)
```

## Arguments

- object:

  A fitted `mlsurv_model` returned by
  [`fit_rsf()`](https://elbadisyimad.com/survalis/reference/fit_rsf.md).

- newdata:

  A `data.frame` of new observations for prediction. Must include the
  same predictor variables used at training time.

- times:

  Optional numeric vector of time points at which to return survival
  probabilities. If `NULL`, predictions are returned at the model's
  internal time grid.

- ...:

  Additional arguments forwarded to
  [`predict.rfsrc`](https://www.randomforestsrc.org//reference/predict.rfsrc.html).

## Value

A base `data.frame` of survival probabilities with one row per
observation in `newdata` and columns named `t={time}` (character),
containing numeric values in \[0, 1\].

## Details

If `times` is provided, the function aligns predictions by selecting,
for each requested time, the closest available time from the RSF
prediction object's `time.interest` grid (nearest-neighbor matching).

## See also

[`fit_rsf()`](https://elbadisyimad.com/survalis/reference/fit_rsf.md),
[`predict.rfsrc`](https://www.randomforestsrc.org//reference/predict.rfsrc.html)

## Examples

``` r
mod_rsf <- fit_rsf(Surv(time, status) ~ age + celltype + karno,
                   data = veteran, ntree = 200)

times <- c(100, 200, 300)
pred_probs <- predict_rsf(mod_rsf, newdata = veteran[1:5, ], times = times)
print(round(pred_probs, 3))
#>   t=100 t=200 t=300
#> 1 0.634 0.313 0.146
#> 2 0.819 0.481 0.282
#> 3 0.672 0.416 0.206
#> 4 0.734 0.448 0.239
#> 5 0.804 0.429 0.244
```
