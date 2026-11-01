# Cross-Validate a Survival Learner (fold-mapped with `fmapn`)

Runs k-fold cross-validation for any pair of `fit_fun`/`pred_fun` that
follow the package's learner contracts, and returns tidy per-fold metric
values. Fold iteration is handled by
[`functionals::fmapn()`](https://rdrr.io/pkg/functionals/man/fmapn.html)
with optional parallel execution and a progress bar.

## Usage

``` r
cv_survlearner(
  formula,
  data,
  fit_fun,
  pred_fun,
  times,
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  verbose = FALSE,
  ncores = 1,
  pb = interactive(),
  ...
)
```

## Arguments

- formula:

  A survival formula `Surv(time, status) ~ predictors`.

- data:

  A data frame containing all variables in `formula`.

- fit_fun:

  Function with signature `fit_fun(formula, data, ...)` that returns an
  `mlsurv_model`.

- pred_fun:

  Function with signature `pred_fun(object, newdata, times, ...)`
  returning a survival-probability matrix/data frame (columns named
  `"t=<time>"`).

- times:

  Numeric vector of evaluation times (passed to `pred_fun`).

- metrics:

  Character vector of metrics to compute. Supported: `"cindex"`,
  `"brier"` (single time), `"ibs"`, `"iae"`, `"ise"`, `"ece"` (single
  time).

- folds:

  Integer; number of folds (default `5`).

- seed:

  Integer random seed for reproducibility (default `123`).

- verbose:

  Logical; print row-dropping due to missingness (default `FALSE`).

- ncores:

  Integer; number of CPU cores for `fmapn` parallel mapping (default
  `1`). Set `> 1` to enable parallelism.

- pb:

  Logical; show a progress bar during fold mapping (default
  [`interactive()`](https://rdrr.io/r/base/interactive.html)).

- ...:

  Additional arguments forwarded to `fit_fun`.

## Value

A `data.frame` with columns: `splits` (rsample split object), `id`,
`fold`, `metric`, and `value`.

## Details

The routine:

1.  Validates `Surv(...)` on the LHS and warns against using `.` in
    formulas.

2.  Drops rows with missing values in any variables referenced by
    `formula`.

3.  Supports `Surv(time, status == k)` by recoding the status to 0/1.

4.  Builds stratified v-folds on the status indicator (rsample).

5.  For each fold: fits on the analysis set, predicts on the assessment
    set, and computes metrics.

Fold iteration is performed via
[`functionals::fmapn()`](https://rdrr.io/pkg/functionals/man/fmapn.html),
which preserves per-fold identifiers (`id`, `fold`) and returns a list
ready for base R's [`rbind()`](https://rdrr.io/r/base/cbind.html).

## See also

[`fmapn`](https://rdrr.io/pkg/functionals/man/fmapn.html)

## Examples

``` r
# \donttest{
cv_results <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno,
  data = veteran,
  fit_fun = fit_coxph,
  pred_fun = predict_coxph,
  times = c(80, 160),
  metrics = c("cindex", "ibs"),
  folds = 2,
  seed = 1
)
# }
```
