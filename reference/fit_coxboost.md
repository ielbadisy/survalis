# Fit a CoxBoost Survival Model (mlsurv_model-compatible)

Fits a likelihood-based boosting model for time-to-event outcomes using
CoxBoost (Binder et al. 2009). Returns an `mlsurv_model` object
compatible with the `survalis` pipeline.

## Usage

``` r
fit_coxboost(
  formula,
  data,
  stepno = 100,
  penalty = NULL,
  criterion = "pscore",
  stepsize.factor = 1
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing variables in `formula`.

- stepno:

  Integer number of boosting steps (default `100`).

- penalty:

  Numeric penalty applied at each boosting step (default
  `9 * number of events`, following `CoxBoost`'s own default).

- criterion:

  Selection criterion for each boosting step: one of `"pscore"`
  (default), `"score"`, `"hpscore"`, `"hscore"`.

- stepsize.factor:

  Numeric step-size factor (default `1`).

## Value

An object of class `mlsurv_model`:

- model:

  Fitted `CoxBoost` model.

- learner:

  `"coxboost"`.

- formula, data:

  Original inputs.

- time, status:

  Names of the survival time and event variables.

- terms:

  Preserved `terms` for consistent prediction design matrices.

## Details

Design contract: returns a named list with engine metadata and preserves
`terms(formula)` so prediction uses the exact same encoding as training.

## Engine

Uses **CoxBoost**, a componentwise boosting algorithm for the Cox
partial likelihood that updates one covariate per boosting step
(penalized score statistic selection), unlike gradient-boosted trees.

## See also

[`predict_coxboost()`](https://elbadisyimad.com/survalis/reference/predict_coxboost.md),
[`tune_coxboost()`](https://elbadisyimad.com/survalis/reference/tune_coxboost.md)

## Examples

``` r
# \donttest{
  mod_cb <- fit_coxboost(
    Surv(time, status) ~ age + karno + celltype,
    data = veteran,
    stepno = 50
  )
  # }
```
