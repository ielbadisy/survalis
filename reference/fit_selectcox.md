# Fit a Predictor-Selection Cox Model (pec::selectCox, mlsurv_model-compatible)

Fits a Cox proportional hazards model with automated predictor selection
using pec's `selectCox()`, returning an `mlsurv_model` object compatible
with the `survalis` evaluation and cross‑validation pipeline.

## Usage

``` r
fit_selectcox(formula, data, rule = "aic")
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`. The
  left‑hand side must be a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object from the
  survival package.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- rule:

  Selection rule passed to
  [`pec::selectCox()`](https://rdrr.io/pkg/pec/man/selectCox.html),
  e.g., `"aic"` (default) or `"p"`.

## Value

An object of class `mlsurv_model`, a named list with elements:

- model:

  The fitted
  [`pec::selectCox`](https://rdrr.io/pkg/pec/man/selectCox.html) model.

- learner:

  Character scalar identifying the learner (`"selectcox"`).

- engine:

  Character scalar naming the engine (`"selectCox"`).

- formula:

  The original survival formula.

- data:

  The training dataset (or a minimal subset needed for prediction).

- rule:

  The selection rule used.

## Details

This wrapper standardizes the return object to the `mlsurv_model`
contract so downstream prediction and evaluation behave consistently
across learners.

## Engine

Uses **pec::selectCox**. Selection can be driven by Akaike's Information
Criterion (`rule = "aic"`) or by p‑value thresholds (`rule = "p"`), as
implemented by pec.

## References

Mogensen UB, Ishwaran H, Gerds TA (2012). Evaluating Random Forests for
Survival Analysis using pec. Gerds TA, et al. pec: Prediction Error
Curves for Survival Models.

## See also

[`predict_selectcox()`](https://elbadisyimad.com/survalis/reference/predict_selectcox.md),
[`tune_selectcox()`](https://elbadisyimad.com/survalis/reference/tune_selectcox.md),
[`pec::selectCox()`](https://rdrr.io/pkg/pec/man/selectCox.html)

## Examples

``` r
mod <- fit_selectcox(Surv(time, status) ~ age + celltype + karno, data = veteran)
```
