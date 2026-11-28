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

groupcv(group, v = 5, seed = NULL)

bootstrap(
  times = 25,
  seed = NULL,
  estimator = c("oob", ".632", ".632+", "optimism"),
  nperm = 20
)
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
  the model formula", whatever its name (for example `cens` in
  `Surv(time, cens) ~ .`). Any other name must be a column of the data;
  otherwise a warning is given and the resampling is not stratified.

- seed:

  Optional integer seed. Set it (and keep it equal across models) for
  paired resampling.

- prop:

  Proportion of rows assigned to the training set in `holdout()`
  (default `0.75`).

- group:

  Column name whose distinct values define the groups held out together
  in `groupcv()`.

- times:

  Integer number of bootstrap resamples (default `25`).

- estimator:

  Bootstrap estimate that
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
  reports as the `mean` of its summary: `"oob"` (out-of-bag, default),
  `".632"`, `".632+"` or `"optimism"`. Any choice other than `"oob"`
  costs one extra fit on all the data and a prediction of all the rows
  by each bootstrap fit.

- nperm:

  Number of random permutations for the no-information score of the
  metrics other than the Brier scores (default `20`).

## Value

An object of class `survalis_resampling`.

## Details

- `cv()` – k-fold cross-validation, optionally repeated, stratified on
  the event indicator by default.

- `holdout()` – a single stratified train/test split.

- `groupcv()` – k-fold cross-validation where whole groups (e.g.
  centres, subjects with repeated rows) are held out together.

- `bootstrap()` – resampling with replacement; the out-of-bag rows form
  the assessment set.
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
  reports the out-of-bag estimate or, with `estimator`, the .632, .632+
  or optimism-corrected estimate (see Details).

Bootstrap estimators. Let \\\overline{err}\\ be the apparent score (the
learner fitted and scored on all the data), \\Err^{(1)}\\ the out-of-bag
score (mean over the resamples of the score of fit \\b\\ on the rows
left out of resample \\b\\) and \\\gamma\\ the no-information score, the
score of the full-data predictions paired with outcomes at random,
averaged over all \\n^2\\ pairs (Efron and Tibshirani 1997, eq. 27). For
the Brier score and the integrated Brier score \\\gamma\\ is computed
exactly at each time (the score is quadratic in the prediction, with
IPCW weights that follow the outcome, as in Gerds and Schumacher 2007);
for the other metrics it is the mean score over `nperm` random
permutations of the prediction rows, an unbiased estimate of the
\\n^2\\-pair mean. Then

- `".632"`: \\0.368\\\overline{err} + 0.632\\Err^{(1)}\\.

- `".632+"`: \\(1 - w)\\\overline{err} + w\\Err^{(1)}\\ with \\w = 0.632
  / (1 - 0.368 R)\\ and relative overfitting rate \\R = (Err^{(1)} -
  \overline{err}) / (\gamma - \overline{err})\\, set to 0 when
  \\Err^{(1)} \le \overline{err}\\ or \\\gamma \le \overline{err}\\ and
  to 1 when \\Err^{(1)} \ge \gamma\\ (Efron and Tibshirani 1997, eqs.
  28-32; the censored-data version of Gerds and Schumacher 2007, Section
  4). With \\R = 1\\ the estimate is \\Err^{(1)}\\.

- `"optimism"`: \\\overline{err}\\ plus the mean over resamples of the
  score of fit \\b\\ on all the data minus its score on its own resample
  (Harrell's optimism bootstrap).

Each estimator is applied to each metric separately (per time for
`"brier@<time>"`, to the integrated score for `"ibs"`). For a metric
where higher is better (`"cindex"`, `"auc"`) the formulas are applied to
its negative, so the inequalities are reversed. The same bootstrap
resamples serve all four estimators, and
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
returns all of them in its `bootstrap` element.
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) and
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
score the out-of-bag rows only.

References: Efron B, Tibshirani R (1997). Improvements on
cross-validation: the .632+ bootstrap method. JASA 92, 548-560. Gerds
TA, Schumacher M (2007). Efron-type measures of prediction error for
survival analysis. Biometrics 63, 1283-1287.

## Examples

``` r
cv(5)
#> <survalis_resampling> 5-fold CV, strata = status
cv(v = 10, repeats = 3, seed = 1)
#> <survalis_resampling> 10-fold CV x 3 repeats, strata = status, seed = 1
holdout(prop = 0.8)
#> <survalis_resampling> holdout (80% train), strata = status
bootstrap(times = 50)
#> <survalis_resampling> bootstrap (50 resamples, oob estimate)
bootstrap(times = 50, estimator = ".632+")
#> <survalis_resampling> bootstrap (50 resamples, .632+ estimate)
```
