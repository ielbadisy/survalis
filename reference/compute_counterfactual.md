# Compute individual counterfactual changes to increase survival

For a single individual (one-row `newdata`), propose feature changes
that maximize survival probability at a target time, subject to optional
per-feature change costs and bounds inferred from the training data in
`model$data`.

## Usage

``` r
compute_counterfactual(
  model,
  newdata,
  times,
  target_time,
  features_to_change = NULL,
  grid.size = 100,
  max.change = NULL,
  cost_penalty = 0.01
)
```

## Arguments

- model:

  A fitted survival model (e.g., an `mlsurv_model`) that was trained
  with a data frame stored in `model$data`, and (ideally) names of the
  outcome variables in `model$time` and `model$status`.

- newdata:

  A data frame with exactly one row representing the individual.

- times:

  Numeric vector of time points used for prediction. Required unless the
  model's predict function can infer times; used together with
  `target_time`.

- target_time:

  Numeric scalar time at which to optimize survival. If missing and
  `times` is provided, the median of `times` is used. If `times` is
  missing but `target_time` is provided, `times` is set to
  `target_time`.

- features_to_change:

  Optional character vector of feature names allowed to change. Defaults
  to all predictors (non-outcome columns).

- grid.size:

  Integer grid size for numeric features (default 100).

- max.change:

  Optional named list of numeric bounds for per-feature absolute change,
  e.g., `list(age = 5)`.

- cost_penalty:

  Numeric penalty weight applied to magnitude of change (default 0.01).

## Value

A `data.frame` with one row per feature considered and columns:
`feature`, `original_value`, `suggested_value`, `survival_gain`,
`change_cost`, `penalized_gain`.

## Details

For each candidate feature, the function sweeps over plausible values
(numeric grid between observed min/max; all other levels for
categorical), predicts survival at `target_time`, and reports the best
penalized gain relative to the original value. Survival predictions are
obtained via a corresponding `predict_*` function inferred from
`model$learner` (e.g., `predict_coxph` for `learner = "coxph"`).

## Examples

``` r
df <- veteran
df$A <- df$trt
mod <- fit_coxph(survival::Surv(time, status) ~ A + age + karno, data = df)

cf <- compute_counterfactual(
  model = mod,
  newdata = df[1, , drop = FALSE],
  times = c(50, 100, 150),
  target_time = 100,
  features_to_change = c("A", "age", "karno"),
  grid.size = 10,
  cost_penalty = 0.02
)
head(cf)
#>   feature original_value suggested_value survival_gain change_cost
#> 1   karno             60         59.4444       -0.0069      0.5556
#> 2     age             69         70.5556        0.0022      1.5556
#> 3       A              1               1        0.0000      0.0000
#>   penalized_gain
#> 1         -0.018
#> 2         -0.029
#> 3          0.000
```
