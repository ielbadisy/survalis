# Performance of a survival learner

Scores survival predictions in one of three ways:

- Resampled:

  `evaluate(fit)` or `evaluate(formula, data, model)` refits the learner
  on each training split of `resampling` and scores it on the held-out
  split (internal validation).

- New data:

  `evaluate(fit, newdata = test)` scores a fitted learner on data it was
  not trained on (external or temporal validation).

- Supplied predictions:

  `evaluate(S, y = Surv(time, status))` scores an `n x length(times)`
  matrix of survival probabilities produced elsewhere, for example by
  another package.

## Usage

``` r
evaluate(object, ...)

# S3 method for class 'survalis_fit'
evaluate(
  object,
  newdata = NULL,
  times = NULL,
  resampling = cv(),
  metrics = c("cindex", "ibs"),
  ncores = 1,
  ...
)

# S3 method for class 'formula'
evaluate(
  object,
  data,
  model,
  spec = list(),
  times = NULL,
  resampling = cv(),
  metrics = c("cindex", "ibs"),
  ncores = 1,
  ...
)

# S3 method for class 'matrix'
evaluate(object, y, times = NULL, metrics = c("cindex", "ibs"), ...)
```

## Arguments

- object:

  A `survalis_fit` from
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md), a
  survival `formula`, or a numeric matrix of survival probabilities
  (rows = subjects, columns = times).

- ...:

  Unused.

- newdata:

  When `object` is a `survalis_fit`: optional data frame with the
  outcome columns. If given, the fit is scored on it once and
  `resampling` is ignored.

- times:

  Evaluation grid. `NULL` picks a grid from the data quantiles and
  records it; for a matrix, `NULL` reads the `"times"` attribute or
  column names of the form `"t=<time>"`.

- resampling:

  A resampling spec from
  [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  or
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  (default
  [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)).

- metrics:

  Character vector from `list_metrics()$metric`. `"brier"` and `"ece"`
  require a single `times` value.

- ncores:

  Cores for the split loop (default `1`).

- data:

  When `object` is a formula: the data frame.

- model:

  When `object` is a formula: a learner id
  (`list_survlearners()$learner`) or a custom learner (see
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md)).

- spec:

  When `object` is a formula: engine arguments (named list).

- y:

  When `object` is a matrix: a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  object with one entry per row of `object`.

## Value

A `survalis_eval` object: `per_split` (long data frame), `summary` (mean
/ sd / se / 95% interval per metric), `times`, `resampling` (`NULL` when
scored once), `metrics`, `model`.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
         spec = list(num.trees = 100))
evaluate(m, times = c(90, 180, 365), resampling = cv(3))
#> <survalis_eval> model = ranger
#> <survalis_resampling> 3-fold CV, strata = status
#>   grid: 3 times in [90, 365]
#>   cindex  0.6410  (sd 0.0790, 95% [0.5520, 0.7300], n = 3)
#>   ibs     0.1790  (sd 0.0180, 95% [0.1590, 0.1990], n = 3)

# train / test split
idx <- seq_len(100)
m_tr <- fit(Surv(time, status) ~ age + karno, veteran[idx, ], "coxph")
evaluate(m_tr, newdata = veteran[-idx, ], times = c(90, 180))
#> <survalis_eval> model = coxph
#>   scored once on 37 subjects (no resampling)
#>   grid: 2 times in [90, 180]
#>   cindex  0.7540
#>   ibs     0.1760

# predictions from anywhere
S <- predict(m_tr, veteran[-idx, ], times = c(90, 180))
evaluate(S, y = Surv(veteran$time[-idx], veteran$status[-idx]))
#> <survalis_eval> model = supplied predictions
#>   scored once on 37 subjects (no resampling)
#>   grid: 2 times in [90, 180]
#>   cindex  0.7540
#>   ibs     0.1760
```
