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
