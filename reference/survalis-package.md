# survalis: Interpretable Survival Machine Learning

survalis provides a unified interface for fitting, predicting,
evaluating, tuning, comparing, interpreting and standardising survival
machine-learning models across 24 learners.

## The seven verbs

The recommended interface, consistent with the funcml package:

- [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md):

  fit a learner by id, returning a `survalis_fit`.

- [`predict()`](https://rdrr.io/r/stats/predict.html):

  typed predictions on a recorded time grid (`"survival"`, `"risk"`,
  `"chf"`, `"hazard"`, `"rmst"`, `"quantile"`, `"median"`).

- [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md):

  resampled performance for a fit or a formula.

- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md):

  grid search, returning the refit best model.

- [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md):

  paired multi-model leaderboard on shared folds.

- [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md):

  one entry point to the explanation methods.

- [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md):

  standardised binary-treatment contrast by G-computation.

Resampling is described by
[`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
and
[`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md).
Discover components with
[`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md),
[`list_metrics()`](https://elbadisyimad.com/survalis/reference/list_metrics.md)
and
[`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md).
Custom learners plug into every verb as
`model = list(fit = , predict = )`; see
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md).

## See also

Useful links:

- <https://elbadisyimad.com/survalis/>

- <https://github.com/ielbadisy/survalis>

- <https://CRAN.R-project.org/package=survalis>

- Report bugs at <https://github.com/ielbadisy/survalis/issues>

## Author

**Maintainer**: Imad El Badisy <elbadisyimad@gmail.com>
