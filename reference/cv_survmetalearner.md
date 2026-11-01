# Cross‑Validate a Stacked Survival Meta‑Learner

Performs v‑fold cross‑validation for the NNLS stacking meta‑learner over
a fixed set of base learners and their predictions. On each fold,
stacking weights are learned on the analysis set and evaluated on the
assessment set.

## Usage

``` r
cv_survmetalearner(
  formula,
  data,
  times,
  base_models,
  base_preds,
  folds = 5,
  metrics = c("cindex", "ibs"),
  seed = 123,
  verbose = TRUE
)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A data frame containing all variables referenced in `formula`.

- times:

  Numeric vector of evaluation times for stacking and scoring.

- base_models:

  A named list of fitted base learner models (used to predict on
  assessment folds and for the final refit on full data).

- base_preds:

  A named list of *training‑set* prediction matrices (rows align with
  `data`, columns align with `times`); names must match `base_models`.

- folds:

  Integer; number of CV folds (default `5`).

- metrics:

  Character vector of metrics to compute (default `c("cindex","ibs")`).
  Supported: `"cindex"`, `"brier"` (single time), `"ibs"`, `"iae"`,
  `"ise"`, `"ece"` (single time).

- seed:

  Integer random seed (default `123`).

- verbose:

  Logical; if `TRUE`, prints fold progress (default `TRUE`).

## Value

An object of class `"cv_survmetalearner_result"` with components:

- model:

  Final `"survmetalearner"` fit on all data.

- cv_results:

  Per‑fold metric values (tibble).

- summary:

  Fold‑aggregated mean and sd by metric (tibble).

- folds, metrics:

  CV settings.

## Details

For each fold: (1) subset `base_preds` to the analysis indices; (2)
learn time‑specific NNLS weights with
[`fit_survmetalearner`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md);

\(3\) predict stacked survival on the assessment set with
[`predict_survmetalearner`](https://elbadisyimad.com/survalis/reference/predict_survmetalearner.md);
(4) compute requested metrics. After CV, a final meta‑learner is fit on
the full data.

## See also

[`fit_survmetalearner`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md),
[`predict_survmetalearner`](https://elbadisyimad.com/survalis/reference/predict_survmetalearner.md),
[`plot_survmetalearner_weights`](https://elbadisyimad.com/survalis/reference/plot_survmetalearner_weights.md)

## Examples

``` r
# \donttest{
form <- Surv(time, status) ~ age + karno + trt
times <- c(80, 160)
mod_cox <- fit_coxph(form, data = veteran)
mod_rpart <- fit_rpart(form, data = veteran)
base_models <- list(coxph = mod_cox, rpart = mod_rpart)
base_preds  <- list(
  coxph = predict_coxph(mod_cox, veteran, times),
  rpart = predict_rpart(mod_rpart, veteran, times)
)

cv_res <- cv_survmetalearner(
  formula = form, data = veteran, times = times,
  base_models = base_models, base_preds = base_preds,
  folds = 2, metrics = c("cindex", "ibs"), seed = 1, verbose = FALSE
)

cv_res$summary
#>   metric   mean        sd
#> 1 cindex 0.7215 0.0516188
#> 2    ibs 0.1790 0.0000000
plot_survmetalearner_weights(cv_res$model)

# }
```
