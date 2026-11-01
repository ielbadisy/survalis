# Changelog

## survalis 1.7.0

- The granular functions (`fit_*()`, `predict_*()`, `tune_*()`,
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md),
  [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md),
  `benchmark*()`, `compute_*()`,
  [`estimate_rmst()`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md),
  …) are no longer scheduled for removal. They stay exported as the
  low-level layer under the verbs and no longer print the one-time
  “granular survalis function” note. `?survalis-deprecated` is now
  `?survalis-lowlevel`.

- [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) no
  longer accepts engine arguments through `...`; pass them in
  `spec = list(...)`. Extra arguments now give an error instead of a
  warning.

- The `stpm2` learner is removed: `fit_stpm2()` and `predict_stpm2()`
  are gone and ‘rstpm2’ is no longer a dependency. Use
  [`fit_rpsurv()`](https://elbadisyimad.com/survalis/reference/fit_rpsurv.md)
  / `fit(model = "rpsurv")`, which fits the same Royston-Parmar model on
  the hazard, odds, or probit scale and adds nonlinear effects (`nle`).

- The `survsvm` learner now has its own engine and no longer depends on
  ‘survivalsvm’. Its “quadprog” path called
  [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html),
  which could stop with “error code 1 from Lapack routine ‘dsyevr’”
  depending on the BLAS/LAPACK library. The new engine fits a kernel
  survival support vector regression on log time (squared loss for
  events, one-sided squared loss for censored times) and solves it
  exactly by active-set Newton steps with one Cholesky factorization
  each, so no quadratic programming solver is needed; it is 10 to 50
  times faster. Survival curves now come from the Kaplan-Meier estimate
  of cross-fitted residuals, `S(t | x) = S_e(log t - f(x))`, which is
  calibrated, instead of `exp(-t / f(x))`, which assumed an exponential
  shape and could give negative times. `dist = "lognormal"` is new;
  `"exp"` and `"weibull"` remain available. The clinical kernel
  (`"add_kernel"`) treats a factor as one predictor and accepts constant
  columns, and `"poly_kernel"` is fixed. Only `type = "regression"` is
  supported (the ranking types give no survival curve); `opt.meth` and
  `diff.meth` are ignored. The default `gamma.mu` is now 0.1 and the
  [`tune_survsvm()`](https://elbadisyimad.com/survalis/reference/tune_survsvm.md)
  grid is `c(0.01, 0.1, 1)`, matched to the new loss scale. The fitted
  object’s `engine` is now `"survalis"`.

- The `bnnsurv` learner (bagged k-nearest-neighbour Kaplan-Meier) now
  has its own engine and no longer depends on ‘bnnSurvival’, which is
  unmaintained. ‘bnnSurvival’ predicted through
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  and forking a session whose multithreaded BLAS (e.g. MKL) is already
  running could hang
  [`predict_bnnsurv()`](https://elbadisyimad.com/survalis/reference/predict_bnnsurv.md).
  The new engine runs in the current process and also: one-hot encodes
  factors (they were integer codes inside the Mahalanobis distance) and
  matches new data to the training levels; uses the pseudo-inverse
  covariance, so constant or collinear columns in a bootstrap sample are
  handled; draws feature subsets from all predictors (a factor counts as
  one); splits weight equally between neighbours tied at the k-th
  distance, so results no longer depend on row order; and reads
  predictions at arbitrary times off the step function instead of
  interpolating linearly. With numeric predictors and default settings
  it gives the same curves as ‘bnnSurvival’ for the same seed, about ten
  times faster. `metric = "euclidean"` (standardized columns) is new.
  The fitted object’s `engine` is now `"survalis"`.

- survalis now declares `R (>= 4.4.0)` and `fastgbm (>= 0.7.0)`. It
  claimed R 4.1, but ‘fastgbm’ up to 0.6.1 needed R 4.5.0 and ‘rms’
  (imported through ‘pec’) needs R 4.4.0, so installation failed on
  older R.

- The `survdnn` learner is removed, and with it the ‘survdnn’ and
  ‘torch’ dependencies. ‘torch’ needs a separate LibTorch download that
  often failed to install. Use the `densemlp` learner
  (`fit(..., model = "densemlp")`), which fits deep survival models with
  a native C++ backend. `fit_survdnn()`, `predict_survdnn()`, and
  `tune_survdnn()` are removed.

- survalis no longer defines a `plot.formula()` method. It replaced
  [`graphics::plot.formula()`](https://rdrr.io/r/graphics/plot.formula.html)
  for every package as soon as survalis was loaded (“Registered S3
  method overwritten”). Use
  [`plot_survcurve()`](https://elbadisyimad.com/survalis/reference/plot_survcurve.md)
  for a Kaplan-Meier curve.

## survalis 1.6.4

Bug fixes found while benchmarking every learner on 31 survival
datasets.

- [`fit_ranger()`](https://elbadisyimad.com/survalis/reference/fit_ranger.md)
  drops predictors that are constant in the data it is fitted on. With
  `respect.unordered.factors = "order"`, ranger stopped with “invalid
  type (NULL) for variable ‘smry\$strata’” on a factor with a single
  level, which happens in tuning folds of data with many rare binary
  predictors (`metabric`).
- [`fit_glmnet()`](https://elbadisyimad.com/survalis/reference/fit_glmnet.md)
  and
  [`predict_glmnet()`](https://elbadisyimad.com/survalis/reference/predict_glmnet.md)
  accept a single predictor. glmnet stops with “x should be a matrix
  with 2 or more columns”, so the predictor matrix is padded with a
  constant zero column.
- [`fit_aalen()`](https://elbadisyimad.com/survalis/reference/fit_aalen.md)
  standardizes numeric covariates and
  [`predict_aalen()`](https://elbadisyimad.com/survalis/reference/predict_aalen.md)
  applies the same scaling.
  [`timereg::aalen()`](https://rdrr.io/pkg/timereg/man/aalen.html)
  silently returned an all-zero fit when covariates were on very
  different scales, so the predicted survival was 1 for every subject.
- [`predict_flexsurvreg()`](https://elbadisyimad.com/survalis/reference/predict_flexsurvreg.md),
  [`predict_bnnsurv()`](https://elbadisyimad.com/survalis/reference/predict_bnnsurv.md),
  and
  [`predict_cforest()`](https://elbadisyimad.com/survalis/reference/predict_cforest.md)
  accept a single evaluation time. A one-time result is now a one-column
  matrix (it was a one-row matrix, or a flat table for `flexsurvreg`),
  so the single-time metrics `brier` and `ece` can be computed for these
  learners.
- [`fit_survmetalearner()`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md)
  can be called with `formula` and `data` only, as
  [`benchmark()`](https://elbadisyimad.com/survalis/reference/benchmark.md)
  and
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  do. It then fits the base learners (`learners`, default `coxph`,
  `glmnet`, `rsf`), builds out-of-fold predictions, and learns the
  weights from them.
  [`predict_survmetalearner()`](https://elbadisyimad.com/survalis/reference/predict_survmetalearner.md)
  interpolates the weights when the requested times differ from the
  fitted grid.
- [`fit_xgboost()`](https://elbadisyimad.com/survalis/reference/fit_xgboost.md)
  accepts `max_depth` and `eta`.
  [`tune_xgboost()`](https://elbadisyimad.com/survalis/reference/tune_xgboost.md)
  passed them, so every configuration failed silently and
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) then
  stopped with an unrelated “Unknown columns: ibs”.
- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) stops
  with a clear message when every configuration fails during
  cross-validation.
- [`tune_survsvm()`](https://elbadisyimad.com/survalis/reference/tune_survsvm.md)
  has a default `param_grid`, and `tune_survdnn()` accepts its own
  default: an argument check rejected the default because
  [`missing()`](https://rdrr.io/r/base/missing.html) is always true for
  it.
- Errors raised inside parallel workers (`ncores > 1`) in
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
  and
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  are now reported.
  [`functionals::fmapn()`](https://rdrr.io/pkg/functionals/man/fmapn.html)
  returns them as `try-error` objects, which were dropped, so a failure
  surfaced as `object 'value' not found`.

## survalis 1.6.3

- Finished the `data.table` -\> `basetable` migration (1.6.1 landed the
  shared `dt-utils.R` choke point; this closes out the rest). All
  remaining
  [`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)/`as.data.table()`/`rbindlist()`/`melt()`/
  `dcast()`/`setorderv()`/`setcolorder()`/`copy()` calls and bare
  `[, list(...), by = ...]`/`:=`/`.SD`/`.N` grouped expressions across
  `list_functions.R`, `plot-survcurve.R`, `compute_shap.R`,
  `compute_shap_mean.R`, `compute_varimp.R`, `compute_pdp.R`,
  `compute_calibration.R`, `fit_survmetalearner.R`, `fit_survdnn.R`,
  `fit_bnnsurv.R`, `fit_coxboost.R`, `fit_xgboost.R`, `fit_glmnet.R`,
  `fit_survsvm.R`, `evaluation.R`
  ([`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md)),
  `survmat-helpers.R`
  ([`plot_survmat()`](https://elbadisyimad.com/survalis/reference/plot_survmat.md)),
  and `benchmark_default_survlearners.R`
  ([`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md)/[`summarize_benchmark_results()`](https://elbadisyimad.com/survalis/reference/summarize_benchmark_results.md)/
  [`best_survlearner()`](https://elbadisyimad.com/survalis/reference/best_survlearner.md))
  are replaced with base R and
  [`basetable::aggregate()`](https://rdrr.io/pkg/basetable/man/aggregate.html)/`count()`/`tolong()`/`towide()`/`orderrows()`/
  `applyby()`. `data.table` dropped from `Imports` entirely.
  [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md)
  and
  [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md)
  etc. now return plain `data.frame`s instead of `data.table`s (no
  downstream code relied on data.table-specific semantics).
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html) grid
  construction is replaced with base
  [`expand.grid()`](https://rdrr.io/r/base/expand.grid.html). All
  existing tests pass; one test asserting a `"data.table"` return class
  was updated to `"data.frame"` to match the new contract.

## survalis 1.6.2

- Fixed a latent correctness bug in the shared `.pmap_rbind_dt()`/
  `.map_rbind_dt()` helpers (introduced in 1.6.1’s basetable migration):
  [`basetable::rbindfill()`](https://rdrr.io/pkg/basetable/man/rbindfill.html)
  silently coerces list-columns (e.g. a per-config vector of
  hidden-layer sizes) to `NA` instead of preserving them, unlike
  [`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html).
  No current learner grid combined through these helpers carries a
  list-valued hyperparameter, so this was not an active regression, but
  it was a landmine for any future one. Both helpers now route through a
  new internal `.rbind_fill_dt()` (base R
  [`rbind()`](https://rdrr.io/r/base/cbind.html) after padding to a
  common column set), which preserves list-columns and still handles
  rows with differing columns (e.g. a failed tuning row that omits
  metric columns).

## survalis 1.6.1

- Migrated the shared `R/dt-utils.R` internal helpers (`.pmap_rbind_dt`,
  `.map_rbind_dt`, `.complete_cases_df`, `.wide_metric_row`,
  `.arrange_by_metric_dt`, `.select_cols`) from `data.table` to
  `basetable` as the row-binding/ordering engine. These helpers are the
  shared choke point behind every tuner’s parameter-grid results table.
  Also fixed a latent single-bracket row-filter in
  [`tune_bnnsurv()`](https://elbadisyimad.com/survalis/reference/tune_bnnsurv.md)
  (`results[!results$failed]`) that only worked under `data.table`’s
  `DT[i]` row-selection semantics; now explicit
  (`results[!results$failed, ]`). `data.table` remains in use elsewhere
  in the package pending further migration.

## survalis 1.6.0

- Bug fix:
  [`survmat_to_rmst()`](https://elbadisyimad.com/survalis/reference/survmat_to_rmst.md)
  (and everything built on it: `predict(type = "rmst")`,
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md))
  previously under-integrated RMST for an off-grid `tau` – it restricted
  the time grid to `times <= tau` and summed whole trapezoidal intervals
  only, silently dropping the fractional last interval instead of
  interpolating into it. For example, on a grid of `0:5`, `tau = 3.5`
  returned the same value as `tau = 3` (wrong). RMST computation is now
  delegated to
  [`tvrmst::rmst_dynamic()`](https://rdrr.io/pkg/tvrmst/man/rmst_dynamic.html),
  with per-subject linear interpolation to `tau` on top of its
  cumulative trapezoidal curve. Values for `tau` that already land
  exactly on a grid point are unchanged; values for off-grid `tau` are
  now correct. `tau` beyond the last grid point still holds at the full
  definite integral (no extrapolation past the observed curve). Added
  `tvrmst` to `Imports`.

## survalis 1.5.0

- New learner: `nbsurv` (conditional naive Bayes survival model, via the
  `nbsurv` package).
  [`fit_nbsurv()`](https://elbadisyimad.com/survalis/reference/fit_nbsurv.md)/[`predict_nbsurv()`](https://elbadisyimad.com/survalis/reference/predict_nbsurv.md)
  added (no tuner). Registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md);
  24 learners total.

## survalis 1.4.0

- New learner: `rpsurv` (Royston-Parmar flexible parametric survival,
  via the `rpsurv` package).
  [`fit_rpsurv()`](https://elbadisyimad.com/survalis/reference/fit_rpsurv.md)/[`predict_rpsurv()`](https://elbadisyimad.com/survalis/reference/predict_rpsurv.md)
  added (no tuner), restricted cubic splines with optional time-varying
  covariate effects (`tve`) and a choice of `scale`
  (hazard/odds/normal). Registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md);
  23 learners total.

## survalis 1.3.0

### Seven-verb interface

- New top-level verbs covering the whole workflow, consistent with the
  `funcml` package:
  - `fit(formula, data, model, spec, seed)` – learner-agnostic fitting;
    returns a `survalis_fit` (a subclass of `mlsurv_model`, so every
    existing helper keeps working on it).
  - `predict(fit, newdata, times, type)` – `type` one of `"survival"`,
    `"risk"`, `"chf"`, `"hazard"` (all returning a validated `survmat`),
    `"rmst"`, `"quantile"`, `"median"` (numeric).
  - `evaluate(fit_or_formula, times, resampling, metrics)` – resampled
    performance; a `survalis_fit` reuses its stored recipe.
  - `tune(formula, data, model, grid, resampling, metric)` – grid
    search, returns the refit best model as a `survalis_fit`.
  - `compare(formula, data, models, specs, resampling, metrics, tune)` –
    paired multi-model leaderboard on shared folds.
  - `interpret(fit, method, newdata, data, times, ...)` – single entry
    point to the `compute_*` explanation methods, with a
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html) that
    dispatches to the paired helper.
  - `estimate(formula, data, model, estimand, treatment, tau)` –
    standardised binary-treatment contrast by G-computation, with a
    subject bootstrap (refit per replicate by default) and an
    identification-assumptions note.

### Supporting infrastructure

- `survmat` prediction contract: an `n x length(times)` numeric matrix
  carrying its grid in the `"times"` attribute (`"t=<time>"` column
  names kept for display and the legacy metric helpers) and the
  predicted quantity in `"quantity"`. Missing predictions and
  non-increasing grids are rejected.
- [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md):
  a recorded, data-dependent evaluation grid used whenever `times` is
  omitted, so runs stay reproducible.
- Resampling constructors
  [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  (with `repeats`),
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md);
  reuse one spec across models for paired comparisons.

### New learners

- `fastgbm` (CRAN): histogram gradient boosting, piecewise-exponential
  objective by default so survival curves are available at arbitrary
  times.
- `densemlp` (CRAN): native-C++ MLP (no `torch`), IPCW-Brier
  discrete-time loss by default, predictions interpolated onto the
  requested times.
- Both registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md);
  22 learners total.

### Deprecations

- The granular functions (`fit_*()`, `predict_*()`, `tune_*()`,
  [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md),
  [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md),
  `benchmark*()`,
  [`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md),
  [`estimate_rmst()`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md),
  …) remain exported and fully supported. They are now marked internal
  in the manual and emit a one-time per-session note. See
  `?\`survalis-deprecated\`\`. They are scheduled to become unexported
  in a future release.

## survalis 1.2.0

- Added
  [`compute_survcluster()`](https://elbadisyimad.com/survalis/reference/compute_survcluster.md)
  and
  [`plot_survcluster()`](https://elbadisyimad.com/survalis/reference/plot_survcluster.md):
  unsupervised clustering of predicted individualized survival curves
  via [`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html)
  (PAM on a weighted curve representation), registered in
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md).
  Added `unsurv` to `Imports` and its paper’s citation (El Badisy 2026,
  `<doi:10.1093/bioadv/vbag218>`) to `Description`.

## survalis 1.1.0

- Added
  [`screen_fdr()`](https://elbadisyimad.com/survalis/reference/screen_fdr.md):
  univariate-Cox variable screening with Benjamini-Hochberg FDR control,
  with a `minscreen` backfill guarantee.
- Added
  [`estimate_rmst()`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md)
  and
  [`plot_rmst()`](https://elbadisyimad.com/survalis/reference/plot_rmst.md):
  standardized RMST estimation via g-computation, with percentile
  bootstrap inference. Dispatch between a marginal estimate and a
  two-arm causal-style contrast is by the `trt_col` argument (`NULL` =
  marginal), matching
  [`plot_survmat()`](https://elbadisyimad.com/survalis/reference/plot_survmat.md)’s
  existing `group` dispatch convention rather than separate function
  names.
- Added
  [`compute_shap_matrix()`](https://elbadisyimad.com/survalis/reference/compute_shap_matrix.md)
  and
  [`plot_shap_beeswarm()`](https://elbadisyimad.com/survalis/reference/plot_shap_beeswarm.md):
  a beeswarm view of SHAP contributions time-integrated over the full
  prediction horizon, built on
  [`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md)’s
  existing per-time Shapley decomposition.
  [`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md)
  gained an explicit `features` argument (defaulting to
  `setdiff(colnames(newdata), c(model$time, model$status))`), replacing
  a hardcoded `c("time", "status", "event")` exclusion, for consistency
  with
  [`compute_interactions()`](https://elbadisyimad.com/survalis/reference/compute_interactions.md)’s
  existing `features` argument.
- Added
  [`roc_survmat()`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  (full time-dependent ROC curve, IPCW-weighted like
  [`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md))
  and
  [`dca_survmat()`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  (decision curve analysis / net benefit for right-censored outcomes,
  using the same `Ĝ(t)` IPCW estimator as
  [`brier()`](https://elbadisyimad.com/survalis/reference/brier.md) and
  [`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)),
  each with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)/[`plot_roc()`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)/[`plot_dca()`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  methods.

## survalis 1.0.0

CRAN release: 2026-08-22

- Added `CoxBoost` as a new survival learner:
  [`fit_coxboost()`](https://elbadisyimad.com/survalis/reference/fit_coxboost.md),
  [`predict_coxboost()`](https://elbadisyimad.com/survalis/reference/predict_coxboost.md),
  and
  [`tune_coxboost()`](https://elbadisyimad.com/survalis/reference/tune_coxboost.md),
  wrapping the `CoxBoost` package’s likelihood-based componentwise
  boosting for the Cox partial likelihood. Registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md)
  and usable directly through
  [`benchmark()`](https://elbadisyimad.com/survalis/reference/benchmark.md).
- Wrapped
  [`tune_glmnet()`](https://elbadisyimad.com/survalis/reference/tune_glmnet.md)’s
  example in `\donttest{}` (per CRAN win-builder pretest feedback: it
  exceeded the 5s/10s example runtime limit on Debian/Windows).

## survalis 0.9.0

- Added a `title` argument to every plotting function
  ([`plot_ale()`](https://elbadisyimad.com/survalis/reference/plot_ale.md),
  [`plot_pdp()`](https://elbadisyimad.com/survalis/reference/plot_pdp.md),
  [`plot_shap()`](https://elbadisyimad.com/survalis/reference/plot_shap.md),
  [`plot_interactions()`](https://elbadisyimad.com/survalis/reference/plot_interactions.md),
  [`plot_calibration()`](https://elbadisyimad.com/survalis/reference/plot_calibration.md),
  [`plot_counterfactual()`](https://elbadisyimad.com/survalis/reference/plot_counterfactual.md),
  [`plot_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_surrogate.md),
  [`plot_tree_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_tree_surrogate.md),
  [`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md),
  [`plot_survmat()`](https://elbadisyimad.com/survalis/reference/plot_survmat.md),
  [`plot_survmetalearner_weights()`](https://elbadisyimad.com/survalis/reference/plot_survmetalearner_weights.md),
  [`plot_benchmark()`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md),
  [`cv_plot()`](https://elbadisyimad.com/survalis/reference/cv_plot.md),
  [`plot_survcurve()`](https://elbadisyimad.com/survalis/reference/plot_survcurve.md)).
  If omitted, the previous automatically generated title is used (fully
  backward compatible); pass `title = NULL` to omit the title entirely,
  for journals that require caption-only figures.

## survalis 0.8.9

- Added a dashed zero-reference line to every plot that can show both
  negative and positive values (previously only
  [`plot_shap()`](https://elbadisyimad.com/survalis/reference/plot_shap.md)’s
  bar chart had one):
  [`plot_ale()`](https://elbadisyimad.com/survalis/reference/plot_ale.md)
  (both `per_time` and `integrated`),
  [`plot_shap()`](https://elbadisyimad.com/survalis/reference/plot_shap.md)’s
  time-dependent curve,
  [`plot_counterfactual()`](https://elbadisyimad.com/survalis/reference/plot_counterfactual.md),
  [`plot_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_surrogate.md),
  and
  [`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md).

## survalis 0.8.8

- [`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md)
  now draws a boxplot of the full per-repetition permutation
  distribution for each feature instead of a single point estimate with
  no uncertainty shown.
  [`compute_varimp()`](https://elbadisyimad.com/survalis/reference/compute_varimp.md)
  retains the per-repetition values (previously discarded after
  computing summary stats) as a `"raw_scores"` attribute on its return
  value, used by
  [`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md);
  falls back to the previous point-plot behavior for objects without
  that attribute (e.g. hand-built summary tables).

## survalis 0.8.7

- Flipped the direction of `plot_interactions(type = "heatmap")`’s
  viridis scale so high interaction values are dark and low values are
  bright, matching the conventional reading of heatmap intensity.

## survalis 0.8.6

- Fixed `plot_interactions(type = "heatmap")`, which was genuinely hard
  to read: the `white -> steelblue` gradient had too little contrast in
  the mid-range, and the diagonal (self-interaction, value 0) blended
  into the low end of the scale. Switched to `scale_fill_viridis_c()`
  (perceptually uniform, high contrast), added an explicit `na.value`
  for any genuinely missing cells, and fixed a squished legend by
  widening the colorbar guide.
- Unified color usage across all plotting functions: grouped plots
  (`color`/`fill` mapped to a variable) previously fell back to
  ggplot2’s default hue palette everywhere except
  [`plot_survcurve()`](https://elbadisyimad.com/survalis/reference/plot_survcurve.md);
  they now all use
  [`scale_color_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md)/[`scale_fill_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md).
  Single-series hardcoded colors (`"steelblue"`, `"skyblue"`,
  `"lightblue"`, `"pink"`, `"tomato"`, plain `"blue"`, and ad-hoc
  green/red hex codes for positive/negative direction) now draw from the
  same shared Dark2-based palette instead of arbitrary named/hex colors.
  This was explicitly deferred in PR
  [\#22](https://github.com/ielbadisy/survalis/issues/22) and is now
  closed out.

## survalis 0.8.5

- Fixed GitHub Actions CI, which had been failing since before this
  development cycle: `man/compute_tree_surrogate.Rd`,
  `man/plot_tree_surrogate.Rd`, and `man/tune_survdnn.Rd` were stale,
  hand-drifted copies using `\donttest{}` (which runs under
  `R CMD check --run-donttest`, as CI does) instead of the `\dontrun{}`
  already present in the roxygen source, causing
  `plot_tree_surrogate(tree_ranger, ...)` to fail on an undefined
  `tree_ranger` object. Regenerated all three from source via
  `devtools::document()`; verified locally with
  `R CMD check --run-donttest` (Status: OK), matching the CI
  configuration.

## survalis 0.8.4

- [`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md)
  now rounds `mean`/`sd`/`se`/`lower`/`upper` to 3 decimals by default
  (`digits = 3` argument), matching
  [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md)’s
  convention; it previously returned full-precision values (e.g.
  `sd = 0.014142136`).

## survalis 0.8.3

- Normalized `base_size` to the shared
  [`theme_survalis()`](https://elbadisyimad.com/survalis/reference/theme_survalis.md)
  default (13) everywhere: `compute_shap.R`, `fit_survmetalearner.R`’s
  weights plot, and
  [`plot_benchmark()`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md)
  previously carried over a hardcoded `base_size = 14` from before the
  theme retrofit (PR
  [\#22](https://github.com/ielbadisy/survalis/issues/22)).

## survalis 0.8.2

- Migrated
  [`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md),
  `compute_shap_mean()`,
  [`compute_calibration()`](https://elbadisyimad.com/survalis/reference/compute_calibration.md),
  and
  [`plot_survmetalearner_weights()`](https://elbadisyimad.com/survalis/reference/plot_survmetalearner_weights.md)
  from dplyr/tidyr to data.table. This was the remaining unqualified
  (non-`dplyr::`-prefixed) dplyr/tidyr usage missed by the original
  file-by-file migration inventory.
- `dplyr`, `tidyr`, `purrr`, and `tibble` are now fully removed from
  `Imports`: the package no longer depends on any tidyverse package for
  its data-manipulation code; `data.table` is the sole engine
  throughout.

## survalis 0.8.1

- Migrated the nested/tuned-CV half of
  `benchmark_default_survlearners.R`
  ([`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md),
  `.nested_surv_*` helpers,
  [`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md),
  [`summarize_benchmark_results()`](https://elbadisyimad.com/survalis/reference/summarize_benchmark_results.md),
  [`best_survlearner()`](https://elbadisyimad.com/survalis/reference/best_survlearner.md))
  from dplyr/tidyr/tibble to data.table. This completes the file-by-file
  dplyr/tidyverse -\> data.table migration across the package. Reuses
  the shared `.score_metrics()` helper in
  `.nested_surv_score_predictions()` instead of a fifth duplicated copy
  of the per-metric scoring
  [`switch()`](https://rdrr.io/r/base/switch.html) block.

## survalis 0.8.0

- Migrated
  [`cv_survmetalearner()`](https://elbadisyimad.com/survalis/reference/cv_survmetalearner.md)
  from dplyr/tidyr/purrr/tibble to data.table, reusing the shared
  `.score_metrics()` helper instead of a third duplicated copy of the
  per-metric scoring [`switch()`](https://rdrr.io/r/base/switch.html)
  block.

## survalis 0.7.27

- Migrated
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md),
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md),
  and
  [`list_metrics()`](https://elbadisyimad.com/survalis/reference/list_metrics.md)
  from
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  to
  [`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)
  (cosmetic descriptive tables, no behavior change).

## survalis 0.7.26

- Migrated
  [`plot_survmat()`](https://elbadisyimad.com/survalis/reference/plot_survmat.md)
  from dplyr/tidyr to data.table
  ([`data.table::melt()`](https://rdrr.io/pkg/data.table/man/melt.data.table.html)
  for the wide-to-long reshape, `[.data.table` with `by =` for the
  group/time summaries).

## survalis 0.7.25

- Migrated
  [`compute_varimp()`](https://elbadisyimad.com/survalis/reference/compute_varimp.md)
  from dplyr/purrr/tibble to data.table.

## survalis 0.7.24

- Migrated `tune_survdnn()` from dplyr/tidyr/purrr/tibble to data.table,
  including the list-valued `hidden` architecture parameter (grid
  expansion built manually over unique levels per parameter, since
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html)
  doesn’t support list-valued columns).

## survalis 0.7.23

- Migrated
  [`tune_survsvm()`](https://elbadisyimad.com/survalis/reference/tune_survsvm.md)
  from dplyr/tidyr/purrr/tibble to data.table (grid expansion via
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html)); the
  CRAN-motivated refit-fallback behavior (PR
  [\#14](https://github.com/ielbadisy/survalis/issues/14)) is unchanged
  and still verified passing.

## survalis 0.7.22

- Migrated
  [`tune_glmnet()`](https://elbadisyimad.com/survalis/reference/tune_glmnet.md)
  from dplyr/tidyr/purrr/tibble to data.table, including its
  [`tidyr::crossing()`](https://tidyr.tidyverse.org/reference/expand.html)
  grid-expansion step (now
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html), a
  sorted-unique cross join with matching semantics).

## survalis 0.7.21

- Migrated
  [`tune_selectcox()`](https://elbadisyimad.com/survalis/reference/tune_selectcox.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.20

- Migrated
  [`tune_flexsurvreg()`](https://elbadisyimad.com/survalis/reference/tune_flexsurvreg.md)
  from dplyr/tidyr/purrr/tibble to data.table. Added `.map_rbind_dt()`
  in `R/dt-utils.R`, a data.table replacement for
  [`purrr::map_dfr()`](https://purrr.tidyverse.org/reference/map_dfr.html)
  over a plain vector (as opposed to `.pmap_rbind_dt()` for a parameter
  grid).

## survalis 0.7.19

- Migrated
  [`tune_blackboost()`](https://elbadisyimad.com/survalis/reference/tune_blackboost.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.18

- Migrated
  [`tune_bart()`](https://elbadisyimad.com/survalis/reference/tune_bart.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.17

- Migrated
  [`tune_rsf()`](https://elbadisyimad.com/survalis/reference/tune_rsf.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.16

- Migrated
  [`tune_orsf()`](https://elbadisyimad.com/survalis/reference/tune_orsf.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.15

- Migrated
  [`tune_cforest()`](https://elbadisyimad.com/survalis/reference/tune_cforest.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.14

- Migrated
  [`tune_bnnsurv()`](https://elbadisyimad.com/survalis/reference/tune_bnnsurv.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.13

- Migrated
  [`tune_xgboost()`](https://elbadisyimad.com/survalis/reference/tune_xgboost.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.12

- Migrated
  [`tune_ranger()`](https://elbadisyimad.com/survalis/reference/tune_ranger.md)
  from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.11

- Fixed
  [`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)
  silently dropping learners whose `tune_*()` has been migrated to
  data.table: `tuning_results[1, cols, drop = FALSE]` relies on
  data.frame `[` semantics, but data.table’s `[` does not select columns
  when `cols` is a variable (only a literal character vector triggers
  that); it was returning the raw column-name vector instead of a
  one-row parameter table, causing every fold to error and the learner
  to be dropped with a warning. Added `.select_cols()` in `R/dt-utils.R`
  that selects columns correctly for either a data.frame/tibble or a
  data.table, so the nested-tuning code stays agnostic to migration
  state.

## survalis 0.7.10

- Retrofitted all existing plotting functions to use
  [`theme_survalis()`](https://elbadisyimad.com/survalis/reference/theme_survalis.md)
  instead of ad-hoc `theme_minimal()` calls with inconsistent
  `base_size` values, so package figures share one consistent visual
  style:
  [`plot_ale()`](https://elbadisyimad.com/survalis/reference/plot_ale.md),
  survmat plotting helpers,
  [`plot_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_surrogate.md)/[`plot_shap()`](https://elbadisyimad.com/survalis/reference/plot_shap.md),
  [`cv_plot()`](https://elbadisyimad.com/survalis/reference/cv_plot.md),
  [`plot_pdp()`](https://elbadisyimad.com/survalis/reference/plot_pdp.md),
  [`plot_interactions()`](https://elbadisyimad.com/survalis/reference/plot_interactions.md),
  [`plot_counterfactual()`](https://elbadisyimad.com/survalis/reference/plot_counterfactual.md),
  [`plot_tree_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_tree_surrogate.md),
  [`plot_benchmark()`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md),
  [`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md),
  [`fit_survmetalearner()`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md)’s
  plotting path. Hardcoded per-plot colors are unchanged in this pass
  (follow-up).

## survalis 0.7.9

- Added
  [`plot_survcurve()`](https://elbadisyimad.com/survalis/reference/plot_survcurve.md),
  a survminer-inspired Kaplan-Meier curve with confidence ribbon,
  log-rank p-value annotation, and an aligned number-at-risk table,
  built natively on
  [`theme_survalis()`](https://elbadisyimad.com/survalis/reference/theme_survalis.md)
  /
  [`scale_color_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md)
  (no dependency on the survminer package). The risk table uses
  `patchwork` (added to Suggests).

## survalis 0.7.8

- Added
  [`theme_survalis()`](https://elbadisyimad.com/survalis/reference/theme_survalis.md),
  a shared ggplot2 theme, and
  [`scale_color_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md)/[`scale_fill_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md),
  colorblind-friendly discrete scales (ColorBrewer “Dark2”-based,
  recycled/interpolated beyond 8 levels), for a consistent visual style
  across package figures. Not yet applied to existing plotting functions
  (follow-up).

## survalis 0.7.7

- Migrated
  [`tune_rpart()`](https://elbadisyimad.com/survalis/reference/tune_rpart.md)
  from dplyr/tidyr/purrr/tibble to data.table. Added internal data.table
  helpers (`R/dt-utils.R`) shared by the ongoing file-by-file migration
  of the remaining `tune_*()`/`fit_*()` grid-search code (tracked in the
  project TODO).

## survalis 0.7.6

- Added
  [`benchmark()`](https://elbadisyimad.com/survalis/reference/benchmark.md),
  the single entry point for comparing survival learners: `tune = FALSE`
  (default) dispatches to
  [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md),
  `tune = TRUE` dispatches to
  [`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)
  (nested CV with per-learner tuning).
  [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md)
  and
  [`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)
  remain available directly.

## survalis 0.7.5

- Migrated the core CV/metrics engine
  ([`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md),
  [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md),
  [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md))
  from dplyr/tidyr/tibble/purrr to data.table; `data.table` is now the
  main data-manipulation engine going forward, replacing the
  dplyr/tidyverse approach used previously. This is a first step; the
  remaining tuning/fitting code (`tune_*()`/`fit_*()`) still uses
  dplyr/purrr and will be migrated incrementally.
- [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md)
  and
  [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md)/[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  metric values now round to 3 decimals by default (`digits = 3`
  argument on
  [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md)).
- These functions now return `data.table` objects instead of tibbles
  ([`is.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) still
  holds; code relying on `tibble`-specific behavior should call
  [`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  explicitly).

## survalis 0.7.4

- Added
  [`timeroc_survmat()`](https://elbadisyimad.com/survalis/reference/timeroc_survmat.md),
  a vectorized cumulative/dynamic time-dependent AUC curve over a vector
  of evaluation times (Uno et al. 2007 / Heagerty and Zheng 2005
  estimator), matching `timeROC::timeROC(weighting = "marginal")` to
  ~1e-3.
- Fixed
  [`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)
  to define cases as events strictly before `t_star` (`time < t_star`)
  rather than `time <= t_star`, matching the canonical Uno/timeROC
  definition; this corrects the `"auc"` metric used throughout
  [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md),
  `benchmark_*()`, and
  [`fit_survmetalearner()`](https://elbadisyimad.com/survalis/reference/fit_survmetalearner.md).
- Added `timeROC` to Suggests for numerical validation in tests.

## survalis 0.7.3

- Fixed
  [`plot_pdp()`](https://elbadisyimad.com/survalis/reference/plot_pdp.md)
  per-time facets using free y-axis scales, which made survival
  probability panels visually incomparable across facets; now uses a
  fixed `[0, 1]` scale via `coord_cartesian()`.

## survalis 0.7.2

- Fixed `tune_survsvm(refit_best = TRUE)` crashing when the top-ranked
  grid candidate fails to refit on the full dataset (e.g., a `quadprog`
  QP infeasibility that only manifests on certain BLAS backends, as seen
  on CRAN’s OpenBLAS check machine); it now falls back to the next-best
  candidate and only errors if every candidate fails.

## survalis 0.7.1

CRAN release: 2026-04-23

- Revised DESCRIPTION to include CRAN-compliant, family-level
  methodological references using validated DOI/URL anchors and
  representative (non-exhaustive) wording.
- Reworked examples across evaluation, benchmarking, interpretation, and
  meta-learner topics to remove commented code and provide executable
  minimal examples.
- Replaced unnecessary `\\dontrun{}` wrappers with `\\donttest{}` for
  executable-but-longer workflows and kept fast examples unwrapped where
  feasible.
- Fixed documentation/example mismatches (including variable-importance
  examples) and regenerated Rd files with `devtools::document()`.
- Verified package checks including `run_dont_test = TRUE` for CRAN
  resubmission readiness.

## survalis 0.7.0

- First CRAN submission candidate.
- Added a cumulative/dynamic time-dependent AUC metric to the evaluation
  layer.
- Added
  [`plot_counterfactual()`](https://elbadisyimad.com/survalis/reference/plot_counterfactual.md)
  for visualizing counterfactual recommendations.
- Added
  [`plot_survmat()`](https://elbadisyimad.com/survalis/reference/plot_survmat.md)
  for plotting predicted survival curves, including grouped summaries.
- Hardened documentation and examples for CRAN checks.
- Updated dependency declarations so exported learners and CV/tuning
  helpers work after installation.
- Added clearer runtime checks for `survdnn` when the LibTorch runtime
  is missing.
