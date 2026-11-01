# Permutation variable importance for survival models

Estimates feature importance by measuring the change in a survival
metric after permuting each feature.

## Usage

``` r
compute_varimp(
  model,
  times,
  metric = "ibs",
  n_repetitions = 10,
  seed = NULL,
  subset = NULL,
  importance_type = c("delta", "mean")
)
```

## Arguments

- model:

  An `mlsurv_model`-like object with fields `data`, `formula`, and
  `learner`.

- times:

  Numeric vector of evaluation times.

- metric:

  Character string, e.g. `"ibs"` or `"cindex"`.

- n_repetitions:

  Integer; number of permutations per feature.

- seed:

  Optional integer seed for reproducibility.

- subset:

  Optional row indices or logical vector to subset `model$data`.

- importance_type:

  One of `"delta"` (default) or `"mean"`; see Details.

## Value

A data.frame with columns:

- `feature`: feature name,

- `importance`: importance value (change in metric),

- `importance_05`, `importance_95`: 5th/95th percentile of the
  per-repetition values used to compute `importance`,

- `scaled_importance`: percent-scaled importance (see Details).

The full per-repetition values (one row per feature per repetition,
columns `feature`, `repetition`, `value`, `scaled_value`) are attached
as the `"raw_scores"` attribute, used by
[`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md)
to draw boxplots instead of a single point per feature.

## Details

For each feature, rows are permuted `n_repetitions` times, predictions
are recomputed, and the chosen metric is compared to the baseline
(unpermuted) value. The `scaled_importance` rescales values to sum to
100%.

## Examples

``` r
# \donttest{
mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
imp <- compute_varimp(
  model = mod,
  times = 80,
  metric = "brier",
  n_repetitions = 3,
  seed = 1,
  subset = 40
)
head(imp)
#> # basetable: 6 x 5
#>    feature importance importance_05 importance_95 scaled_importance
#> 1    karno 0.10133333        0.0762        0.1239        100.000000
#> 2 celltype 0.04566667        0.0152        0.0710         45.065789
#> 3      age 0.00200000        0.0011        0.0029          1.973684
#> 4      trt 0.00000000        0.0000        0.0000          0.000000
#> 5 diagtime 0.00000000        0.0000        0.0000          0.000000
#> 6    prior 0.00000000        0.0000        0.0000          0.000000
# }
```
