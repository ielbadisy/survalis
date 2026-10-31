## Update from 1.0.0 to 1.7.0

Main changes (see NEWS.md for the full list):

* A seven-verb interface (`fit()`, `predict()`, `evaluate()`, `tune()`,
  `compare()`, `interpret()`, `estimate()`) on top of the existing learners.
* New learners: `rpsurv`, `nbsurv`, `fastgbm`, `densemlp`.
* Removed dependencies: 'torch' and 'survdnn' (the `survdnn` learner is
  defunct, replaced by `densemlp`) and 'bnnSurvival', which is unmaintained;
  the `bnnsurv` learner now has its own engine. 'bnnSurvival' predicted
  through `parallel::mclapply()`, which could hang when forked from a
  session running a multithreaded BLAS.
* `fit_survsvm()` refits with `opt.meth = "ipop"` when the "quadprog" solver
  fails, which depends on the BLAS/LAPACK library.
* survalis no longer registers a `plot.formula()` method that overrode
  `graphics::plot.formula()`.
* `Depends: R (>= 4.4.0)` and `fastgbm (>= 0.7.0)`: the earlier R (>= 4.1)
  was wrong ('rms', imported through 'pec', needs R 4.4.0).

This version needs 'fastgbm' 0.7.0, submitted to CRAN just before it.

## R CMD check results

0 errors | 0 warnings | 1 note

* The note is the "Imports includes 35 non-default packages" INFO: survalis
  wraps learners from many packages behind one interface.

## Test environments

* local Ubuntu 24.04, R 4.5.1 (MKL BLAS), `R CMD check --as-cran`
