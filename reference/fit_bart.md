# Fit a Bayesian Additive Regression Trees (BART) Survival Model

Fits a right‐censored survival model using BART's
[`surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html) and returns
an `mlsurv_model` compatible with the survalis pipeline.

## Usage

``` r
fit_bart(formula, data, K = 3, ...)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`. Only
  right‐censored responses are supported.

- data:

  A `data.frame` containing variables referenced in `formula`.

- K:

  Integer; number of internal time grid intervals used by the BART
  survival engine (default `3`). Larger values yield a finer internal
  evaluation grid.

- ...:

  Further args to BART::surv.bart

## Value

An object of class `mlsurv_model` with elements:

- model:

  Fitted
  [`BART::surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html)
  object.

- learner:

  `"bart"`.

- formula, data:

  Original inputs.

- eval_times:

  Engine's internal time grid used for prediction.

## Details

The response must be of class `Surv(..., type = "right")`. The fitted
object stores the engine's internal evaluation times in `$eval_times`
(from `bart_fit$times`) for downstream prediction alignment.

## See also

[`predict_bart`](https://elbadisyimad.com/survalis/reference/predict_bart.md),
[`surv.bart`](https://rdrr.io/pkg/BART/man/surv.bart.html)

## Examples

``` r
ex_data <- veteran[1:40, c("time", "status", "age", "karno", "celltype")]
mod_bart <- fit_bart(
  Surv(time, status) ~ age + karno + celltype,
  data = ex_data,
  K = 1,
  ntree = 5,
  ndpost = 20,
  nskip = 5,
  mc.cores = 1,
  seed = 42
)
```
