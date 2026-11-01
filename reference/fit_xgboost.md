# Fit an XGBoost Survival Model (mlsurv_model-compatible)

Fits a gradient-boosted tree model for time-to-event outcomes using
xgboost. Supports `"survival:aft"` (default) and `"survival:cox"`.
Returns an `mlsurv_model` object compatible with the `survalis`
pipeline.

## Usage

``` r
fit_xgboost(
  formula,
  data,
  booster = "gbtree",
  objective = "survival:aft",
  aft_loss_distribution = "extreme",
  aft_loss_distribution_scale = 1,
  nrounds = 100,
  max_depth = 6,
  eta = 0.3
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing variables in `formula`.

- booster:

  XGBoost booster type (default `"gbtree"`).

- objective:

  One of `"survival:aft"` (default) or `"survival:cox"`.

- aft_loss_distribution:

  AFT error distribution: `"extreme"` (default), `"logistic"`, or
  `"normal"`.

- aft_loss_distribution_scale:

  Positive numeric scale for the AFT loss (default `1`).

- nrounds:

  Integer number of boosting iterations (default `100`).

- max_depth:

  Maximum tree depth (default `6`).

- eta:

  Learning rate (default `0.3`).

## Value

An object of class `mlsurv_model`:

- model:

  Fitted `xgboost` model.

- learner:

  `"xgboost"`.

- formula, data:

  Original inputs.

- time, status:

  Names of the survival time and event variables.

- objective:

  Objective used.

- dist, scale:

  AFT distribution and scale (populated even if Cox is used).

- terms:

  Preserved `terms` for consistent prediction design matrices.

## Details

Design contract: returns a named list with engine metadata and preserves
`terms(formula)` so prediction uses the exact same encoding as training.

## Engine

Uses **xgboost**. For `"survival:aft"`, interval labels are set via
`label_lower_bound`/`label_upper_bound` with right-censoring represented
by `Inf` on the upper bound. For `"survival:cox"`, labels follow the Cox
objective convention.

## See also

[`predict_xgboost()`](https://elbadisyimad.com/survalis/reference/predict_xgboost.md),
[`tune_xgboost()`](https://elbadisyimad.com/survalis/reference/tune_xgboost.md)

## Examples

``` r
# \donttest{
  mod_xgb <- fit_xgboost(
    Surv(time, status) ~ age + karno + celltype,
    data = veteran
  )
  # }
```
