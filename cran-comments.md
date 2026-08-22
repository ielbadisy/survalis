## Resubmission of an archived package

`survalis` 0.7.1 was archived from CRAN on 2026-08-22 (issues not corrected
in time: `tune_survsvm(refit_best = TRUE)` could crash on `quadprog::solve.QP`
QP infeasibility under R-devel/OpenBLAS on CRAN's check machine). That fix
was already made in 0.7.2 but was never actually submitted to CRAN before
the archival deadline; it is included in this 1.0.0 submission (see
`tune_survsvm()`'s candidate-fallback logic in `R/fit_survsvm.R`).

This submission also fixes a win-builder pretest failure from an earlier
1.0.0 submission attempt: `tune_glmnet()`'s example exceeded CRAN's example
runtime limit (10.64s on Windows, 5.68s on Debian; limits are 10s/5s
respectively). Fixed by wrapping the example in `\donttest{}`, matching the
existing convention already used for other tuning-function examples in this
package (`tune_orsf()`, etc.).

## R CMD check results

0 errors | 0 warnings | see below for notes

## Reverse dependencies

None known.
