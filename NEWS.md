# survalis 1.3.0

## Seven-verb interface

* New top-level verbs covering the whole workflow, consistent with the
  `funcml` package:
  * `fit(formula, data, model, spec, seed)` -- learner-agnostic fitting;
    returns a `survalis_fit` (a subclass of `mlsurv_model`, so every existing
    helper keeps working on it).
  * `predict(fit, newdata, times, type)` -- `type` one of `"survival"`,
    `"risk"`, `"chf"`, `"hazard"` (all returning a validated `survmat`),
    `"rmst"`, `"quantile"`, `"median"` (numeric).
  * `evaluate(fit_or_formula, times, resampling, metrics)` -- resampled
    performance; a `survalis_fit` reuses its stored recipe.
  * `tune(formula, data, model, grid, resampling, metric)` -- grid search,
    returns the refit best model as a `survalis_fit`.
  * `compare(formula, data, models, specs, resampling, metrics, tune)` --
    paired multi-model leaderboard on shared folds.
  * `interpret(fit, method, newdata, data, times, ...)` -- single entry point
    to the `compute_*` explanation methods, with a `plot()` that dispatches to
    the paired helper.
  * `estimate(formula, data, model, estimand, treatment, tau)` -- standardised
    binary-treatment contrast by G-computation, with a subject bootstrap
    (refit per replicate by default) and an identification-assumptions note.

## Supporting infrastructure

* `survmat` prediction contract: an `n x length(times)` numeric matrix
  carrying its grid in the `"times"` attribute (`"t=<time>"` column names kept
  for display and the legacy metric helpers) and the predicted quantity in
  `"quantity"`. Missing predictions and non-increasing grids are rejected.
* `default_times()`: a recorded, data-dependent evaluation grid used whenever
  `times` is omitted, so runs stay reproducible.
* Resampling constructors `cv()` (with `repeats`), `holdout()`, `group_cv()`,
  `bootstrap()`; reuse one spec across models for paired comparisons.

## New learners

* `fastgbm` (CRAN): histogram gradient boosting, piecewise-exponential
  objective by default so survival curves are available at arbitrary times.
* `densemlp` (CRAN): native-C++ MLP (no `torch`), IPCW-Brier discrete-time
  loss by default, predictions interpolated onto the requested times.
* Both registered in `list_survlearners()`; 22 learners total.

## Deprecations

* The granular functions (`fit_*()`, `predict_*()`, `tune_*()`,
  `cv_survlearner()`, `score_survmodel()`, `benchmark*()`,
  `summarise_benchmark()`, `estimate_rmst()`, ...) remain exported and fully
  supported. They are now marked internal in the manual and emit a one-time
  per-session note. See `?\`survalis-deprecated\``. They are scheduled to
  become unexported in a future release.

# survalis 1.2.0

* Added `compute_survcluster()` and `plot_survcluster()`: unsupervised
  clustering of predicted individualized survival curves via
  `unsurv::unsurv()` (PAM on a weighted curve representation), registered
  in `list_interpretability_methods()`. Added `unsurv` to `Imports` and
  its paper's citation (El Badisy 2026, `<doi:10.1093/bioadv/vbag218>`)
  to `Description`.

# survalis 1.1.0

* Added `screen_fdr()`: univariate-Cox variable screening with
  Benjamini-Hochberg FDR control, with a `minscreen` backfill guarantee.
* Added `estimate_rmst()` and `plot_rmst()`: standardized RMST
  estimation via g-computation, with percentile bootstrap inference. Dispatch
  between a marginal estimate and a two-arm causal-style contrast is by the
  `trt_col` argument (`NULL` = marginal), matching `plot_survmat()`'s
  existing `group` dispatch convention rather than separate function names.
* Added `compute_shap_matrix()` and `plot_shap_beeswarm()`: a beeswarm view
  of SHAP contributions time-integrated over the full prediction horizon,
  built on `compute_shap()`'s existing per-time Shapley decomposition.
  `compute_shap()` gained an explicit `features` argument (defaulting to
  `setdiff(colnames(newdata), c(model$time, model$status))`), replacing a
  hardcoded `c("time", "status", "event")` exclusion, for consistency with
  `compute_interactions()`'s existing `features` argument.
* Added `roc_survmat()` (full time-dependent ROC curve, IPCW-weighted like
  `auc_survmat()`) and `dca_survmat()` (decision curve analysis / net benefit
  for right-censored outcomes, using the same `Ĝ(t)` IPCW estimator as
  `brier()` and `auc_survmat()`), each with `plot()`/`plot_roc()`/`plot_dca()`
  methods.

# survalis 1.0.0

* Added `CoxBoost` as a new survival learner: `fit_coxboost()`,
  `predict_coxboost()`, and `tune_coxboost()`, wrapping the `CoxBoost`
  package's likelihood-based componentwise boosting for the Cox partial
  likelihood. Registered in `list_survlearners()` and usable directly
  through `benchmark()`.
* Wrapped `tune_glmnet()`'s example in `\donttest{}` (per CRAN win-builder
  pretest feedback: it exceeded the 5s/10s example runtime limit on
  Debian/Windows).

# survalis 0.9.0

* Added a `title` argument to every plotting function (`plot_ale()`,
  `plot_pdp()`, `plot_shap()`, `plot_interactions()`, `plot_calibration()`,
  `plot_counterfactual()`, `plot_surrogate()`, `plot_tree_surrogate()`,
  `plot_varimp()`, `plot_survmat()`, `plot_survmetalearner_weights()`,
  `plot_benchmark()`, `cv_plot()`, `plot_survcurve()`). If omitted, the
  previous automatically generated title is used (fully backward
  compatible); pass `title = NULL` to omit the title entirely, for
  journals that require caption-only figures.

# survalis 0.8.9

* Added a dashed zero-reference line to every plot that can show both
  negative and positive values (previously only `plot_shap()`'s bar chart
  had one): `plot_ale()` (both `per_time` and `integrated`), `plot_shap()`'s
  time-dependent curve, `plot_counterfactual()`, `plot_surrogate()`, and
  `plot_varimp()`.

# survalis 0.8.8

* `plot_varimp()` now draws a boxplot of the full per-repetition permutation
  distribution for each feature instead of a single point estimate with no
  uncertainty shown. `compute_varimp()` retains the per-repetition values
  (previously discarded after computing summary stats) as a `"raw_scores"`
  attribute on its return value, used by `plot_varimp()`; falls back to the
  previous point-plot behavior for objects without that attribute (e.g.
  hand-built summary tables).

# survalis 0.8.7

* Flipped the direction of `plot_interactions(type = "heatmap")`'s viridis
  scale so high interaction values are dark and low values are bright,
  matching the conventional reading of heatmap intensity.

# survalis 0.8.6

* Fixed `plot_interactions(type = "heatmap")`, which was genuinely hard to
  read: the `white -> steelblue` gradient had too little contrast in the
  mid-range, and the diagonal (self-interaction, value 0) blended into the
  low end of the scale. Switched to `scale_fill_viridis_c()` (perceptually
  uniform, high contrast), added an explicit `na.value` for any genuinely
  missing cells, and fixed a squished legend by widening the colorbar guide.
* Unified color usage across all plotting functions: grouped plots
  (`color`/`fill` mapped to a variable) previously fell back to ggplot2's
  default hue palette everywhere except `plot_survcurve()`; they now all use
  `scale_color_survalis()`/`scale_fill_survalis()`. Single-series hardcoded
  colors (`"steelblue"`, `"skyblue"`, `"lightblue"`, `"pink"`, `"tomato"`,
  plain `"blue"`, and ad-hoc green/red hex codes for positive/negative
  direction) now draw from the same shared Dark2-based palette instead of
  arbitrary named/hex colors. This was explicitly deferred in PR #22 and is
  now closed out.

# survalis 0.8.5

* Fixed GitHub Actions CI, which had been failing since before this
  development cycle: `man/compute_tree_surrogate.Rd`, `man/plot_tree_surrogate.Rd`,
  and `man/tune_survdnn.Rd` were stale, hand-drifted copies using
  `\donttest{}` (which runs under `R CMD check --run-donttest`, as CI does)
  instead of the `\dontrun{}` already present in the roxygen source,
  causing `plot_tree_surrogate(tree_ranger, ...)` to fail on an undefined
  `tree_ranger` object. Regenerated all three from source via
  `devtools::document()`; verified locally with
  `R CMD check --run-donttest` (Status: OK), matching the CI configuration.

# survalis 0.8.4

* `summarise_benchmark()` now rounds `mean`/`sd`/`se`/`lower`/`upper` to 3
  decimals by default (`digits = 3` argument), matching `cv_summary()`'s
  convention; it previously returned full-precision values (e.g.
  `sd = 0.014142136`).

# survalis 0.8.3

* Normalized `base_size` to the shared `theme_survalis()` default (13)
  everywhere: `compute_shap.R`, `fit_survmetalearner.R`'s weights plot, and
  `plot_benchmark()` previously carried over a hardcoded `base_size = 14`
  from before the theme retrofit (PR #22).

# survalis 0.8.2

* Migrated `compute_shap()`, `compute_shap_mean()`, `compute_calibration()`,
  and `plot_survmetalearner_weights()` from dplyr/tidyr to data.table. This
  was the remaining unqualified (non-`dplyr::`-prefixed) dplyr/tidyr usage
  missed by the original file-by-file migration inventory.
* `dplyr`, `tidyr`, `purrr`, and `tibble` are now fully removed from
  `Imports`: the package no longer depends on any tidyverse package for
  its data-manipulation code; `data.table` is the sole engine throughout.

# survalis 0.8.1

* Migrated the nested/tuned-CV half of `benchmark_default_survlearners.R`
  (`benchmark_default_survlearners()`, `.nested_surv_*` helpers,
  `summarise_benchmark()`, `summarize_benchmark_results()`,
  `best_survlearner()`) from dplyr/tidyr/tibble to data.table. This
  completes the file-by-file dplyr/tidyverse -> data.table migration across
  the package. Reuses the shared `.score_metrics()` helper in
  `.nested_surv_score_predictions()` instead of a fifth duplicated copy of
  the per-metric scoring `switch()` block.

# survalis 0.8.0

* Migrated `cv_survmetalearner()` from dplyr/tidyr/purrr/tibble to
  data.table, reusing the shared `.score_metrics()` helper instead of a
  third duplicated copy of the per-metric scoring `switch()` block.

# survalis 0.7.27

* Migrated `list_survlearners()`, `list_interpretability_methods()`, and
  `list_metrics()` from `tibble::tibble()` to `data.table::data.table()`
  (cosmetic descriptive tables, no behavior change).

# survalis 0.7.26

* Migrated `plot_survmat()` from dplyr/tidyr to data.table
  (`data.table::melt()` for the wide-to-long reshape, `[.data.table` with
  `by =` for the group/time summaries).

# survalis 0.7.25

* Migrated `compute_varimp()` from dplyr/purrr/tibble to data.table.

# survalis 0.7.24

* Migrated `tune_survdnn()` from dplyr/tidyr/purrr/tibble to data.table,
  including the list-valued `hidden` architecture parameter (grid expansion
  built manually over unique levels per parameter, since `data.table::CJ()`
  doesn't support list-valued columns).

# survalis 0.7.23

* Migrated `tune_survsvm()` from dplyr/tidyr/purrr/tibble to data.table
  (grid expansion via `data.table::CJ()`); the CRAN-motivated refit-fallback
  behavior (PR #14) is unchanged and still verified passing.

# survalis 0.7.22

* Migrated `tune_glmnet()` from dplyr/tidyr/purrr/tibble to data.table,
  including its `tidyr::crossing()` grid-expansion step (now
  `data.table::CJ()`, a sorted-unique cross join with matching semantics).

# survalis 0.7.21

* Migrated `tune_selectcox()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.20

* Migrated `tune_flexsurvreg()` from dplyr/tidyr/purrr/tibble to data.table.
  Added `.map_rbind_dt()` in `R/dt-utils.R`, a data.table replacement for
  `purrr::map_dfr()` over a plain vector (as opposed to `.pmap_rbind_dt()`
  for a parameter grid).

# survalis 0.7.19

* Migrated `tune_blackboost()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.18

* Migrated `tune_bart()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.17

* Migrated `tune_rsf()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.16

* Migrated `tune_orsf()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.15

* Migrated `tune_cforest()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.14

* Migrated `tune_bnnsurv()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.13

* Migrated `tune_xgboost()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.12

* Migrated `tune_ranger()` from dplyr/tidyr/purrr/tibble to data.table.

# survalis 0.7.11

* Fixed `benchmark_tuned_survlearners()` silently dropping learners whose
  `tune_*()` has been migrated to data.table: `tuning_results[1, cols, drop = FALSE]`
  relies on data.frame `[` semantics, but data.table's `[` does not select
  columns when `cols` is a variable (only a literal character vector
  triggers that); it was returning the raw column-name vector instead of
  a one-row parameter table, causing every fold to error and the learner
  to be dropped with a warning. Added `.select_cols()` in `R/dt-utils.R`
  that selects columns correctly for either a data.frame/tibble or a
  data.table, so the nested-tuning code stays agnostic to migration state.

# survalis 0.7.10

* Retrofitted all existing plotting functions to use `theme_survalis()`
  instead of ad-hoc `theme_minimal()` calls with inconsistent `base_size`
  values, so package figures share one consistent visual style: `plot_ale()`,
  survmat plotting helpers, `plot_surrogate()`/`plot_shap()`, `cv_plot()`,
  `plot_pdp()`, `plot_interactions()`, `plot_counterfactual()`,
  `plot_tree_surrogate()`, `plot_benchmark()`, `plot_varimp()`,
  `fit_survmetalearner()`'s plotting path. Hardcoded per-plot colors are
  unchanged in this pass (follow-up).

# survalis 0.7.9

* Added `plot_survcurve()`, a survminer-inspired Kaplan-Meier curve with
  confidence ribbon, log-rank p-value annotation, and an aligned
  number-at-risk table, built natively on `theme_survalis()` /
  `scale_color_survalis()` (no dependency on the survminer package). The
  risk table uses `patchwork` (added to Suggests).

# survalis 0.7.8

* Added `theme_survalis()`, a shared ggplot2 theme, and
  `scale_color_survalis()`/`scale_fill_survalis()`, colorblind-friendly
  discrete scales (ColorBrewer "Dark2"-based, recycled/interpolated beyond
  8 levels), for a consistent visual style across package figures. Not yet
  applied to existing plotting functions (follow-up).

# survalis 0.7.7

* Migrated `tune_rpart()` from dplyr/tidyr/purrr/tibble to data.table.
  Added internal data.table helpers (`R/dt-utils.R`) shared by the ongoing
  file-by-file migration of the remaining `tune_*()`/`fit_*()` grid-search
  code (tracked in the project TODO).

# survalis 0.7.6

* Added `benchmark()`, the single entry point for comparing survival
  learners: `tune = FALSE` (default) dispatches to
  `benchmark_default_survlearners()`, `tune = TRUE` dispatches to
  `benchmark_tuned_survlearners()` (nested CV with per-learner tuning).
  `benchmark_default_survlearners()` and `benchmark_tuned_survlearners()`
  remain available directly.

# survalis 0.7.5

* Migrated the core CV/metrics engine (`cv_survlearner()`, `cv_summary()`,
  `score_survmodel()`) from dplyr/tidyr/tibble/purrr to data.table; `data.table`
  is now the main data-manipulation engine going forward, replacing the
  dplyr/tidyverse approach used previously. This is a first step; the
  remaining tuning/fitting code (`tune_*()`/`fit_*()`) still uses dplyr/purrr
  and will be migrated incrementally.
* `cv_summary()` and `score_survmodel()`/`cv_survlearner()` metric values now
  round to 3 decimals by default (`digits = 3` argument on `cv_summary()`).
* These functions now return `data.table` objects instead of tibbles
  (`is.data.frame()` still holds; code relying on `tibble`-specific behavior
  should call `tibble::as_tibble()` explicitly).

# survalis 0.7.4

* Added `timeroc_survmat()`, a vectorized cumulative/dynamic time-dependent
  AUC curve over a vector of evaluation times (Uno et al. 2007 / Heagerty
  and Zheng 2005 estimator), matching `timeROC::timeROC(weighting = "marginal")`
  to ~1e-3.
* Fixed `auc_survmat()` to define cases as events strictly before `t_star`
  (`time < t_star`) rather than `time <= t_star`, matching the canonical
  Uno/timeROC definition; this corrects the `"auc"` metric used throughout
  `score_survmodel()`, `benchmark_*()`, and `fit_survmetalearner()`.
* Added `timeROC` to Suggests for numerical validation in tests.

# survalis 0.7.3

* Fixed `plot_pdp()` per-time facets using free y-axis scales, which made
  survival probability panels visually incomparable across facets; now uses
  a fixed `[0, 1]` scale via `coord_cartesian()`.

# survalis 0.7.2

* Fixed `tune_survsvm(refit_best = TRUE)` crashing when the top-ranked grid
  candidate fails to refit on the full dataset (e.g., a `quadprog` QP
  infeasibility that only manifests on certain BLAS backends, as seen on
  CRAN's OpenBLAS check machine); it now falls back to the next-best
  candidate and only errors if every candidate fails.

# survalis 0.7.1

* Revised DESCRIPTION to include CRAN-compliant, family-level methodological references using validated DOI/URL anchors and representative (non-exhaustive) wording.
* Reworked examples across evaluation, benchmarking, interpretation, and meta-learner topics to remove commented code and provide executable minimal examples.
* Replaced unnecessary `\\dontrun{}` wrappers with `\\donttest{}` for executable-but-longer workflows and kept fast examples unwrapped where feasible.
* Fixed documentation/example mismatches (including variable-importance examples) and regenerated Rd files with `devtools::document()`.
* Verified package checks including `run_dont_test = TRUE` for CRAN resubmission readiness.

# survalis 0.7.0

* First CRAN submission candidate.
* Added a cumulative/dynamic time-dependent AUC metric to the evaluation layer.
* Added `plot_counterfactual()` for visualizing counterfactual recommendations.
* Added `plot_survmat()` for plotting predicted survival curves, including grouped summaries.
* Hardened documentation and examples for CRAN checks.
* Updated dependency declarations so exported learners and CV/tuning helpers work after installation.
* Added clearer runtime checks for `survdnn` when the LibTorch runtime is missing.
