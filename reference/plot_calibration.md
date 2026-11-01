# Plot Calibration Curve for Survival Predictions

Produces a calibration plot comparing mean predicted survival (x-axis)
to observed survival with bootstrap CIs (y-axis) at a single evaluation
time.

## Usage

``` r
plot_calibration(calib_output, smooth = TRUE, title)
```

## Arguments

- calib_output:

  Output list returned by
  [`compute_calibration()`](https://elbadisyimad.com/survalis/reference/compute_calibration.md).

- smooth:

  Logical; if `TRUE`, overlays a LOESS smooth of
  `observed_surv ~ mean_pred_surv`.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A `ggplot2` object showing bin-wise calibration points, bootstrap error
bars, the 45° reference line, and (optionally) a smooth curve.

## Details

Points above the diagonal indicate underprediction (observed survival
higher than predicted), while points below indicate overprediction.

## See also

[`compute_calibration()`](https://elbadisyimad.com/survalis/reference/compute_calibration.md)

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
calib <- compute_calibration(
  model = mod,
  data = veteran,
  time = "time",
  status = "status",
  eval_time = 80,
  n_bins = 4,
  n_boot = 5,
  seed = 1
)
plot_calibration(calib)
#> Warning: span too small.   fewer data values than degrees of freedom.
#> Warning: pseudoinverse used at 0.13942
#> Warning: neighborhood radius 0.46201
#> Warning: reciprocal condition number  0
#> Warning: There are other near singularities as well. 0.12856

```
