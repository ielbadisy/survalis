# Resampled performance of a survival learner

Refits `model` on each training split of `resampling` and scores its
predictions on the held-out split, returning per-split metric values
plus a summary. A `survalis_fit` may be passed directly, in which case
its stored formula, data and `spec` are reused.

## Usage

``` r
evaluate(object, ...)

# S3 method for class 'survalis_fit'
evaluate(
  object,
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
```

## Arguments

- object:

  A `survalis_fit` from
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md), or a
  survival `formula`.

- ...:

  Unused.

- times:

  Evaluation grid; `NULL` uses
  [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md)
  and records it.

- resampling:

  A
  [resampling](https://elbadisyimad.com/survalis/reference/resampling.md)
  spec (default
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
  (`list_survlearners()$learner`).

- spec:

  When `object` is a formula: engine arguments (named list).

## Value

A `survalis_eval` object: `per_split` (long data frame), `summary` (mean
/ sd / se / 95% interval per metric), `times`, `resampling`, `metrics`,
`model`.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
         spec = list(num.trees = 100))
evaluate(m, times = c(90, 180, 365), resampling = cv(3))
#> <survalis_eval> model = ranger
#> <survalis_resampling> 3-fold CV, strata = status
#>   grid: 3 times in [90, 365]
#>   cindex  0.5870  (sd 0.0560, 95% [0.5240, 0.6500], n = 3)
#>   ibs     0.1710  (sd 0.0350, 95% [0.1310, 0.2100], n = 3)
```
