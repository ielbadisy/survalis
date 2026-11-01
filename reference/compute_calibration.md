# Calibration of Survival Predictions at a Single Time Point

Computes a nonparametric calibration curve for a survival model at one
evaluation time by binning predicted survival probabilities and
comparing bin-wise means to Kaplan-Meier-based observed survival, with
bootstrap CIs.

## Usage

``` r
compute_calibration(
  model,
  data,
  time,
  status,
  eval_time,
  n_bins = 10,
  n_boot = 100,
  seed = 123,
  learner_name = NULL
)
```

## Arguments

- model:

  An `mlsurv_model` fitted via a `fit_*()` function; must include a
  valid `learner` so the corresponding `predict_<learner>()` can be
  dispatched, and (ideally) a `$data` field for level alignment.

- data:

  A data frame with predictors and survival outcome columns.

- time:

  Survival time; either a numeric vector of the same length as
  `nrow(data)` or a single string giving the column name in `data`.

- status:

  Event indicator; either a numeric/logical vector or a single string
  giving the column name in `data`. Events should be coded 1, censoring
  0.

- eval_time:

  Single numeric time at which to assess calibration.

- n_bins:

  Integer number of quantile-based bins used to group predictions.

- n_boot:

  Integer number of bootstrap resamples for confidence intervals.

- seed:

  Integer seed for reproducibility of binning/bootstrap.

- learner_name:

  Optional character override for labeling the learner in downstream
  plots (defaults to `model$learner`).

## Value

A list with components:

- calibration_table:

  A data frame with columns `bin`, `mean_pred_surv`, `observed_surv`,
  `lower_ci`, `upper_ci`.

- eval_time:

  The evaluation time used.

- n_bins:

  Number of bins.

- n_boot:

  Number of bootstrap resamples.

- learner:

  The learner label (from `learner_name` or `model$learner`).

## Details

Predicted survival at `eval_time` is obtained from the appropriate
`predict_<learner>()`. Predictions are split into `n_bins` quantile
bins. For each bin, the function reports: mean predicted survival,
observed survival at `eval_time` from a Kaplan-Meier fit on the bin's
rows, and bootstrap percentile (2.5%, 97.5%) CIs on the observed
survival computed by resampling rows with replacement.

## See also

[`plot_calibration()`](https://elbadisyimad.com/survalis/reference/plot_calibration.md)

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
head(calib$calibration_table)
#>   bin mean_pred_surv observed_surv   lower_ci  upper_ci
#> 1   1      0.1425391    0.05714286 0.03003861 0.1470779
#> 2   2      0.4101197    0.47058824 0.26256109 0.5330159
#> 3   3      0.6014329    0.61764706 0.51559252 0.7023308
#> 4   4      0.7655523    0.84926471 0.74432900 0.9122807
```
