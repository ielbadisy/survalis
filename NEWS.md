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
  triggers that) — it was returning the raw column-name vector instead of
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
