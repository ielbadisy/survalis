# Fit a Stacked Survival Meta‑Learner (Time‑Varying NNLS)

Learns nonnegative, time‑specific stacking weights over a set of base
survival learners by solving a nonnegative least squares (NNLS) problem
at each requested time \\t^{\*}\\. For each \\t^{\*}\\, the target is
the indicator \\Y = I(T \> t^{\*})\\; the features are the base
learners' predicted survival probabilities \\S\_\ell(t^{\*} \mid x)\\.
Weights are constrained to be nonnegative and are normalized to sum to 1
per time point.

## Usage

``` r
fit_survmetalearner(
  base_preds = NULL,
  time = NULL,
  status = NULL,
  times = NULL,
  base_models = NULL,
  formula = NULL,
  data = NULL,
  learners = c("coxph", "glmnet", "rsf"),
  folds = 5,
  seed = 123L
)
```

## Arguments

- base_preds:

  A named list of matrices/data frames, one per base learner, each of
  dimension `n x length(times)`, with columns named `"t=<time>"`. When
  `NULL` (the default), the base learners in `learners` are fitted from
  `formula` and `data`, and `base_preds` are their out-of-fold
  predictions.

- time:

  Numeric vector of observed event/censoring times (length `n`). Taken
  from `data` when `base_preds` is `NULL`.

- status:

  Numeric/binary vector of event indicators (1=event, 0=censor) (length
  `n`). Taken from `data` when `base_preds` is `NULL`.

- times:

  Numeric vector of evaluation times at which to learn weights. When
  `NULL` and `base_preds` is `NULL`, a grid of 20 times from
  [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md)
  is used.

- base_models:

  A named list of fitted base learner objects; names must match
  `names(base_preds)` and the learner names used by `predict_*()`.
  Fitted automatically when `base_preds` is `NULL`.

- formula:

  A [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) formula that
  was used for the base models (stored for metadata and downstream
  scoring).

- data:

  The training data frame used for the base models (stored for
  metadata).

- learners:

  Base learners fitted when `base_preds` is `NULL` (default
  `c("coxph", "glmnet", "rsf")`). A learner that fails in cross-fitting
  is dropped with a warning.

- folds:

  Number of cross-fitting folds used to build the out-of-fold
  predictions (default `5`).

- seed:

  Seed for the cross-fitting folds.

## Value

An object of class `c("mlsurv_model","survmetalearner")` with elements:

- weights:

  Matrix `L x T` of nonnegative stacking weights, rows=learners,
  cols=`"t=<time>"`.

- base_models:

  Named list of base learner fits.

- base_preds:

  Named list of base prediction matrices on training data.

- learners:

  Character vector of learner names (from `names(base_preds)`).

- formula, data, time, status:

  Training metadata for scoring/reporting.

- learner:

  The string `"survmetalearner"` (for `predict_*` dispatch).

## Details

For each `t` in `times`, this function fits \$\$\min\_{w \ge 0} \\mid
Y - S w \\\_2^2 \quad \text{s.t.} \\ \sum_j w_j = 1\$\$ where \\Y = I(T
\> t)\\ and \\S\\ is the matrix of base survival probabilities at \\t\\.
The NNLS solution from nnls is renormalized to sum to 1 (if the solution
is all zeros, weights remain `NA` for that time).

## See also

[`predict_survmetalearner`](https://elbadisyimad.com/survalis/reference/predict_survmetalearner.md),
[`plot_survmetalearner_weights`](https://elbadisyimad.com/survalis/reference/plot_survmetalearner_weights.md),
[`cv_survmetalearner`](https://elbadisyimad.com/survalis/reference/cv_survmetalearner.md)

## Examples

``` r
form <- Surv(time, status) ~ age + karno + trt
times <- c(80, 160)
mod_cox <- fit_coxph(form, data = veteran)
mod_rpart <- fit_rpart(form, data = veteran)
base_models <- list(coxph = mod_cox, rpart = mod_rpart)
base_preds <- list(
  coxph = predict_coxph(mod_cox, newdata = veteran, times = times),
  rpart = predict_rpart(mod_rpart, newdata = veteran, times = times)
)
meta_model <- fit_survmetalearner(
  base_preds = base_preds,
  time = veteran$time,
  status = veteran$status,
  times = times,
  base_models = base_models,
  formula = form,
  data = veteran
)
meta_model$weights
#>            t=80 t=160
#> coxph 0.6630347     0
#> rpart 0.3369653     1
```
