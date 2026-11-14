## Resubmission

This resubmission fixes the incoming pre-test NOTE "Base packages in
Suggests/Enhances imported in NAMESPACE": 'stats' and 'utils' moved from
Suggests to Imports.

The NOTE "Possibly misspelled words in DESCRIPTION: El, Badisy" is a false
positive: it is the author's surname in a citation.

## Update from 1.0.0 to 2.4.1

Main changes (see NEWS.md for the full list):

* A seven-verb interface (`fit()`, `predict()`, `evaluate()`, `tune()`,
  `compare()`, `interpret()`, `estimate()`) on top of the existing learners.
  The exported API shrinks from 108 to 16 functions: the granular `fit_*()`,
  `predict_*()`, `compute_*()`, `plot_*()` and metric helpers are now
  internal, mapped to their verb in `vignette("from-1.0")`. survalis has no
  reverse dependencies on CRAN.
* The remaining exported helpers are renamed to glued lowercase, like the
  verbs: `list_survlearners()`, `list_metrics()`,
  `list_interpretability_methods()` and `group_cv()` become `listlearners()`,
  `listmetrics()`, `listmethods()` and `groupcv()`. The underscore names are
  removed, with no wrappers; the old-to-new table is in NEWS.md (2.4.0).
* New learners: `rpsurv`, `nbsurv`, `fastgbm`, `densemlp`. New
  decision-weighted scores (`dwscore()`), Kaplan-Meier, cumulative
  incidence and model-adjusted curves (`kmplot()`), and `interpret()` methods `"roc"`, `"dca"` and
  `"survcluster"`.
* Removed dependencies: 'torch' and 'survdnn' (the `survdnn` learner is
  removed, replaced by `densemlp`) and 'bnnSurvival', which is unmaintained;
  the `bnnsurv` learner now has its own engine. 'bnnSurvival' predicted
  through `parallel::mclapply()`, which could hang when forked from a
  session running a multithreaded BLAS.
* Removed dependency 'survivalsvm': the `survsvm` learner now has its own
  engine with no quadratic programming solver. The 'survivalsvm' path called
  `Matrix::nearPD()`, which failed with a LAPACK error on some BLAS/LAPACK
  libraries.
* Removed dependencies 'rstpm2' (the `stpm2` learner is replaced by
  `rpsurv`) and 'data.table' (replaced by 'basetable'). New dependencies,
  all on CRAN: 'basetable' (>= 1.4.5), 'unsurv', 'tvrmst', 'fastgbm', 'densemlp',
  'rpsurv', 'nbsurv'. Internal table operations (reshaping, grouped
  summaries, row binding, sorting) now use 'basetable'.
* survalis no longer registers a `plot.formula()` method that overrode
  `graphics::plot.formula()`.
* `Depends: R (>= 4.4.0)` and `fastgbm (>= 0.7.0)`: the earlier R (>= 4.1)
  was wrong ('rms', imported through 'pec', needs R 4.4.0). 'fastgbm' 0.7.0
  is now on CRAN.

## R CMD check results

0 errors | 0 warnings | 0 notes

* `R CMD check` reports the INFO "Imports includes 32 non-default packages":
  survalis wraps learners from many packages behind one interface.

## Test environments

* local Ubuntu 24.04, R 4.5.1 (MKL BLAS), `R CMD check --as-cran`
