# Resampling specifications for evaluate(), tune() and compare()

These constructors return a lightweight `survalis_resampling`
description that the verb-level functions turn into concrete train/test
splits. They do not touch the data themselves, so the same specification
can be reused across models to obtain paired comparisons (identical
folds via a shared `seed`).

## Usage

``` r
cv(v = 5, repeats = 1, strata = "status", seed = NULL)

holdout(prop = 0.75, strata = "status", seed = NULL)

group_cv(group, v = 5, seed = NULL)

bootstrap(times = 25, seed = NULL)
```

## Arguments

- v:

  Integer number of folds (default `5`).

- repeats:

  Integer number of times the whole v-fold scheme is repeated with fresh
  random folds (default `1`).

- strata:

  Column name to stratify on, or `NA` to disable stratification. The
  default sentinel `"status"` means "use the event-indicator column of
  the model formula", resolved by the caller.

- seed:

  Optional integer seed. Set it (and keep it equal across models) for
  paired resampling.

- prop:

  Proportion of rows assigned to the training set in `holdout()`
  (default `0.75`).

- group:

  Column name whose distinct values define the groups held out together
  in `group_cv()`.

- times:

  Integer number of bootstrap resamples (default `25`).

## Value

An object of class `survalis_resampling`.

## Details

- `cv()` – k-fold cross-validation, optionally repeated, stratified on
  the event indicator by default.

- `holdout()` – a single stratified train/test split.

- `group_cv()` – k-fold cross-validation where whole groups (e.g.
  centres, subjects with repeated rows) are held out together.

- `bootstrap()` – resampling with replacement; the out-of-bag rows form
  the assessment set (Efron-Tibshirani style, suitable for `.632`
  reporting).

## Examples

``` r
cv(5)
#> <survalis_resampling> 5-fold CV, strata = status
cv(v = 10, repeats = 3, seed = 1)
#> <survalis_resampling> 10-fold CV x 3 repeats, strata = status, seed = 1
holdout(prop = 0.8)
#> <survalis_resampling> holdout (80% train), strata = status
bootstrap(times = 50)
#> <survalis_resampling> bootstrap (50 resamples, OOB assessment)
```
