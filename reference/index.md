# Package index

## Verbs

Fit a learner, predict, evaluate, compare, tune, interpret and estimate
contrasts from one formula and one dataset. Every learner, built in or
your own, goes in `model`.

- [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) : Fit a
  survival learner
- [`predict(`*`<survalis_fit>`*`)`](https://elbadisyimad.com/survalis/reference/predict.survalis_fit.md)
  : Predict from a fitted survival learner
- [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
  : Performance of a survival learner
- [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  : Compare several survival learners on shared resampling folds
- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) : Tune
  a survival learner over a hyperparameter grid
- [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
  : Interpret a fitted survival learner
- [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)
  : Standardised causal-style contrast for a binary treatment
- [`as.data.frame(`*`<survmat>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_fit>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_eval>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_tune>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_compare>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_estimate>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  : Tidy data frames from verb result objects
- [`survalis`](https://elbadisyimad.com/survalis/reference/survalis-package.md)
  [`survalis-package`](https://elbadisyimad.com/survalis/reference/survalis-package.md)
  : survalis: Interpretable Survival Machine Learning

## Resampling

Split schemes used by
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) and
[`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md).

- [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  : Resampling specifications for evaluate(), tune() and compare()

## Scores

Decision-weighted scores over stated risk thresholds.

- [`dwscore()`](https://elbadisyimad.com/survalis/reference/dwscore.md)
  : Decision-weighted scores for survival predictions

## Lookups

Valid ids for `model`, `metrics` and `method`.

- [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md)
  : List the survival learners
- [`list_metrics()`](https://elbadisyimad.com/survalis/reference/list_metrics.md)
  : List the evaluation metrics
- [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md)
  : List the interpretability methods

## Data and re-exports

- [`veteran`](https://elbadisyimad.com/survalis/reference/veteran.md) :
  Veteran's Administration Lung Cancer Trial Data
- [`reexports`](https://elbadisyimad.com/survalis/reference/reexports.md)
  [`Surv`](https://elbadisyimad.com/survalis/reference/reexports.md) :
  Re-export Surv from survival
