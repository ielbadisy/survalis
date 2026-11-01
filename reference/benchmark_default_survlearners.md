# Benchmark Multiple Survival Learners (Cross-Validation Wrapper)

Runs
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
for a set of learner names (e.g., `"ranger"`, `"coxph"`) by dynamically
dispatching `fit_<learner>` and `predict_<learner>` functions. Returns
the row‑bound CV results across all requested learners.

## Usage

``` r
benchmark_default_survlearners(
  formula,
  data,
  learners,
  times,
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ncores = 1,
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
  `c("ranger","coxph","glmnet")`. For each id, `fit_<id>` and
  `predict_<id>` must exist.

- times:

  Numeric vector of evaluation time points for survival predictions.

- metrics:

  Character vector of metrics to compute in CV (e.g., `"cindex"`,
  `"ibs"`, `"iae"`, `"ise"`). The metric semantics are those of
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md).

- folds:

  Integer number of CV folds (default `5`).

- seed:

  Integer random seed for fold generation.

- ncores:

  Integer; number of CPU cores passed to
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  for fold-level parallel mapping (default `1`).

- verbose:

  Logical; if `TRUE`, prints progress per learner.

- suppress_errors:

  Logical; if `TRUE` (default) errors from individual learners are
  caught and reported via
  [`warning()`](https://rdrr.io/r/base/warning.html) and benchmarking
  continues. If `FALSE`, the first error is thrown.

- ...:

  Additional arguments forwarded to each learner's `fit_*()`.

## Value

A data frame of CV results (as returned by
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md))
with an extra column `learner` identifying the source learner.

## Details

Learners are run independently using identical CV splits and scoring
settings. Any learner whose `fit_*()` or `predict_*()` function is
missing will be skipped with a warning. At least one learner must
complete successfully or an error is raised.

## See also

[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md),
[`plot_benchmark()`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md),
[`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md)

## Examples

``` r
# \donttest{
res <- benchmark_default_survlearners(
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
