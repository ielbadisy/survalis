# Compare several survival learners on shared resampling folds

Evaluates each learner in `models` on the *same* set of resampling
splits, so the per-fold values are paired across models. Returns a
leaderboard plus the full per-split table.

## Usage

``` r
compare(
  formula,
  data,
  models,
  specs = list(),
  times = NULL,
  resampling = cv(),
  metrics = c("cindex", "ibs"),
  tune = FALSE,
  ncores = 1,
  ...
)
```

## Arguments

- formula:

  A survival formula.

- data:

  A data frame.

- models:

  Character vector of learner ids.

- specs:

  Named list of engine-argument lists, keyed by model id
  (`specs[["ranger"]] <- list(num.trees = 500)`).

- times:

  Evaluation grid; `NULL` uses
  [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md).

- resampling:

  A
  [resampling](https://elbadisyimad.com/survalis/reference/resampling.md)
  spec (default
  [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md));
  the split set is drawn once and reused for every model.

- metrics:

  Character vector from `list_metrics()$metric`.

- tune:

  Logical. When `TRUE`, each tunable model is tuned once on the full
  data (its default grid) and the winning configuration is then
  evaluated on the shared folds. This is not nested CV, so the tuned
  results are optimistic; a message is emitted. Default `FALSE`.

- ncores:

  Cores for the split loop (default `1`).

- ...:

  Unused.

## Value

A `survalis_compare` object: `leaderboard` (model x metric, mean and
se), `per_split` (long), `times`, `resampling`, `metrics`, `models`,
`tuned`.

## Examples

``` r
cmp <- compare(Surv(time, status) ~ age + karno + celltype, veteran,
               models = c("coxph", "ranger"), times = c(90, 180, 365),
               resampling = cv(3))
cmp
#> <survalis_compare> 2 models
#> <survalis_resampling> 3-fold CV, strata = status
#> 
#>   cindex (higher is better):
#>    * coxph            0.7433  (se 0.0305)
#>      ranger           0.6600  (se 0.0165)
#> 
#>   ibs (lower is better):
#>    * coxph            0.1493  (se 0.0234)
#>      ranger           0.1600  (se 0.0193)
plot(cmp)
```
