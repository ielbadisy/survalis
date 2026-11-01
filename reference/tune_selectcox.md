# Tune SelectCox Rule (Cross-Validation)

Cross‑validates pec's `selectCox()` across one or more selection rules
and selects the best configuration by the primary metric. Optionally
refits the best rule on the full dataset.

## Usage

``` r
tune_selectcox(
  formula,
  data,
  times,
  rules = c("aic", "p"),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ncores = 1,
  refit_best = TRUE,
  ...
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`. The
  left‑hand side must be a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) object from the
  survival package.

- data:

  A `data.frame` containing all variables referenced in `formula`.

- times:

  Numeric vector of evaluation time points (same scale as the training
  survival time). Must be non‑negative and finite.

- rules:

  Character vector of selection rules to compare (e.g., `c("aic","p")`).

- metrics:

  Character vector of metrics to evaluate/optimize (e.g., `"cindex"`,
  `"ibs"`, `"ise"`). The first entry is the primary selection metric.

- folds:

  Number of cross‑validation folds.

- seed:

  Integer random seed for reproducibility.

- ncores:

  Integer number of CPU cores passed to
  [`cv_survlearner`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold evaluation (default `1`).

- refit_best:

  Logical; if `TRUE`, refits the best rule on the full data and returns
  it as an `mlsurv_model` (with tuning results attached).

- ...:

  Additional arguments forwarded to the underlying engine where
  applicable.

## Value

If `refit_best = FALSE`, a `data.frame` (class `"tuned_surv"`) with a
row per rule and metric columns, ordered by the first metric. If
`refit_best = TRUE`, a fitted `mlsurv_model` from
[`fit_selectcox()`](https://elbadisyimad.com/survalis/reference/fit_selectcox.md)
with attribute `"tuning_results"` containing the full results.

## Details

**Evaluation.** Internally calls
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
with
[`fit_selectcox()`](https://elbadisyimad.com/survalis/reference/fit_selectcox.md)/[`predict_selectcox()`](https://elbadisyimad.com/survalis/reference/predict_selectcox.md)
to ensure tuning uses the same code paths as production. Rules typically
include `"aic"` and/or `"p"`.

## See also

[`fit_selectcox()`](https://elbadisyimad.com/survalis/reference/fit_selectcox.md),
[`predict_selectcox()`](https://elbadisyimad.com/survalis/reference/predict_selectcox.md),
[`pec::selectCox()`](https://rdrr.io/pkg/pec/man/selectCox.html)

## Examples

``` r
res_selectcox <- tune_selectcox(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  rules = c("aic", "p"),
  metrics = c("cindex", "ibs", "ise"),
  folds = 3,
  refit_best = FALSE
)

print(res_selectcox)
#> # basetable: 2 x 4
#>   rule cindex   ibs    ise
#> 1  aic  0.725 0.152 14.348
#> 2    p  0.725 0.152 14.348
class(res_selectcox)
#> [1] "tuned_surv" "basetable"  "data.frame"

mod_selectcox <- tune_selectcox(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  rules = c("aic", "p"),
  metrics = c("cindex", "ibs", "ise"),
  folds = 3,
  refit_best = TRUE
)

summary(mod_selectcox)
#> 
#> ── selectcox summary ───────────────────────────────────────────────────────────
#> Formula:
#> Surv(time, status) ~ age + karno + celltype
#> Engine: selectCox
#> Learner: selectcox
#> Data summary:
#> - Observations: 137
#> - Predictors: "age, karno, celltypesmallcell, celltypeadeno, celltypelarge"
#> - Time range: [1, 999]
#> - Event rate: "93.4%"
```
