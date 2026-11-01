# Accumulated Local Effects (ALE) for Survival Models

Computes ALE curves for a numeric (continuous) feature with respect to
survival probabilities at one or more evaluation times. ALE summarizes
the average *local* effect of changing a feature within small intervals,
is robust to correlated features, and is centered to have mean zero.

## Usage

``` r
compute_ale(model, newdata, feature, times, grid.size = 20)
```

## Arguments

- model:

  An `mlsurv_model` created by a `fit_*()` function. Must include a
  valid `learner` so that `predict_<learner>()` can be dispatched.

- newdata:

  Data frame used to compute ALE (typically the training set or a
  representative sample).

- feature:

  Single **numeric/continuous** feature name for which to compute ALE.
  Categorical features are not supported here (use PDP/ICE).

- times:

  Numeric vector of time points at which to evaluate survival
  probabilities.

- grid.size:

  Integer number of quantile cut points used to build the ALE grid
  (default 20). The algorithm uses quantiles of `newdata[[feature]]`.

## Value

A list with:

- ale:

  A data frame with columns `feature_value` and one column per time
  (`"t=<time>"`) containing centered ALE effects.

- integrated:

  If multiple times were provided, a data frame with columns
  `feature_value` and `integrated_ale` equal to the mean of per-time ALE
  effects across `times`; otherwise `NULL`.

## Details

For consecutive quantile bins \\\[z_k, z\_{k+1}\]\\ of the target
feature, ALE integrates the *local* change in the model prediction when
moving the feature from \\z_k\\ to \\z\_{k+1}\\ while holding all other
features at their observed values, and then accumulates these
differences across bins. For survival models, predictions are survival
probabilities at `times`. The returned ALE curves are centered (mean
zero across the grid) per time.

## See also

[`plot_ale()`](https://elbadisyimad.com/survalis/reference/plot_ale.md),
[`compute_pdp()`](https://elbadisyimad.com/survalis/reference/compute_pdp.md)

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
ale_res <- compute_ale(
  model = mod,
  newdata = veteran,
  feature = "karno",
  times = c(80, 160),
  grid.size = 8
)
head(ale_res$ale)
#>   feature_value        t=80       t=160
#> 1            10 -0.40564642 -0.24056434
#> 2            30 -0.27222595 -0.21351483
#> 3            50 -0.06932490 -0.11616708
#> 4            60  0.03751524 -0.03392164
#> 5            70  0.13682692  0.06526931
#> 6            80  0.22367683  0.17252028
```
