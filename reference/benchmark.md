# Benchmark Survival Learners (Meta Function)

The single entry point for comparing multiple survival learners. With
`tune = FALSE` (default), each learner is run with fixed hyperparameters
via
[`benchmark_default_survlearners`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md).
With `tune = TRUE`, each learner is tuned internally per outer fold via
nested cross-validation using
[`benchmark_tuned_survlearners`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md).

## Usage

``` r
benchmark(
  formula,
  data,
  learners,
  times,
  metrics = c("cindex", "ibs"),
  tune = FALSE,
  folds = 5,
  outer_folds = 5,
  inner_folds = 5,
  seed = 123,
  ncores = 1,
  inner_ncores = 1,
  learner_args = list(),
  refit_final = FALSE,
  verbose = FALSE,
  suppress_errors = TRUE,
  ...
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ x1 + x2 + ...`.

- data:

  A data frame containing the variables in `formula`.

- learners:

  Character vector of learner ids (without prefixes), e.g.
  `c("ranger", "coxph", "glmnet")`.

- times:

  Numeric vector of evaluation time points for survival predictions.

- metrics:

  Character vector of metrics to compute/optimize.

- tune:

  Logical; if `FALSE` (default), benchmarks learners with fixed
  hyperparameters (`folds`, `ncores` apply). If `TRUE`, tunes each
  learner via nested CV (`outer_folds`, `inner_folds`, `inner_ncores`,
  `learner_args`, `refit_final` apply).

- folds:

  Integer number of CV folds; used only when `tune = FALSE`.

- outer_folds, inner_folds:

  Integer numbers of outer/inner CV folds; used only when `tune = TRUE`.

- seed:

  Integer random seed.

- ncores:

  Integer CPU cores for fold-level mapping; used only when
  `tune = FALSE`.

- inner_ncores:

  Integer CPU cores for each learner's inner tuning CV; used only when
  `tune = TRUE`.

- learner_args:

  Optional named list of learner-specific tuning arguments; used only
  when `tune = TRUE`. See
  [`benchmark_tuned_survlearners`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md).

- refit_final:

  Logical; if `TRUE` (and `tune = TRUE`), tunes and refits each learner
  on the full dataset in addition to the outer CV.

- verbose:

  Logical; print progress.

- suppress_errors:

  Logical; if `TRUE` (default), per-learner errors are caught and
  reported via [`warning()`](https://rdrr.io/r/base/warning.html)
  instead of stopping the whole benchmark.

- ...:

  Additional arguments forwarded to
  [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md)
  (fit arguments) when `tune = FALSE`, or to
  [`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)
  (common tuning arguments) when `tune = TRUE`.

## Value

With `tune = FALSE`, a data frame of CV results (see
[`benchmark_default_survlearners`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md)).
With `tune = TRUE`, a list of class `"nested_surv_benchmark"` (see
[`benchmark_tuned_survlearners`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)).

## See also

[`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md),
[`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)

## Examples

``` r
# \donttest{
res <- benchmark(
  Surv(time, status) ~ age + karno + trt,
  data = veteran,
  learners = c("coxph", "rpart"),
  times = c(80, 160),
  metrics = c("cindex", "ibs"),
  folds = 2,
  seed = 1
)
head(res)
#>         splits    id fold metric value learner
#> 1 list(trt.... Fold1    1 cindex 0.694   coxph
#> 2 list(trt.... Fold1    1    ibs 0.197   coxph
#> 3 list(trt.... Fold2    2 cindex 0.714   coxph
#> 4 list(trt.... Fold2    2    ibs 0.184   coxph
#> 5 list(trt.... Fold1    1 cindex 0.650   rpart
#> 6 list(trt.... Fold1    1    ibs 0.270   rpart
# }
```
