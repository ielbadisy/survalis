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
