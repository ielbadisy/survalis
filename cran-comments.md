## Resubmission of an archived package

`survalis` 0.7.1 was archived from CRAN on 2026-08-22 (issues not corrected
in time: `tune_survsvm(refit_best = TRUE)` could crash on `quadprog::solve.QP`
QP infeasibility under R-devel/OpenBLAS on CRAN's check machine). That fix
was already made in 0.7.2 but was never actually submitted to CRAN before
the archival deadline; it was included in the 1.0.0 submission (see
`tune_survsvm()`'s candidate-fallback logic in `R/fit_survsvm.R`) and remains
in place.

The 1.0.0 submission also fixed a win-builder pretest failure from an earlier
1.0.0 submission attempt: `tune_glmnet()`'s example exceeded CRAN's example
runtime limit (10.64s on Windows, 5.68s on Debian; limits are 10s/5s
respectively), fixed by wrapping the example in `\donttest{}`, matching the
existing convention already used for other tuning-function examples in this
package (`tune_orsf()`, etc.).

## New in this submission (1.1.0)

Four additive features, all built on `survalis`'s existing survmat contract
and IPCW machinery (no new dependencies):

* `screen_fdr()`: univariate-Cox variable screening with Benjamini-Hochberg
  FDR control.
* `estimate_rmst()` / `plot_rmst()`: standardized RMST estimation via
  g-computation, with a `trt_col` argument dispatching between a marginal
  estimate and a two-arm contrast (percentile bootstrap inference).
* `compute_shap_matrix()` / `plot_shap_beeswarm()`: a beeswarm view of SHAP
  contributions integrated over the full prediction horizon.
* `roc_survmat()` / `dca_survmat()`: time-dependent ROC and decision curve
  analysis for right-censored outcomes, using the package's existing IPCW
  censoring-weight estimator.

All new `@examples` run under `R CMD check`; heavier ones are wrapped in
`\donttest{}` following the existing convention. Full `testthat` coverage
added for each new function.

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the automated CRAN incoming-feasibility message:

```
Days since last update: 3
```

There are no package-code, documentation, example, test, or manual notes.

## Reverse dependencies

None known.
