# survalis 0.7.0

* First CRAN submission candidate.
* Added a cumulative/dynamic time-dependent AUC metric to the evaluation layer.
* Added `plot_counterfactual()` for visualizing counterfactual recommendations.
* Added `plot_survmat()` for plotting predicted survival curves, including grouped summaries.
* Hardened documentation and examples for CRAN checks.
* Updated dependency declarations so exported learners and CV/tuning helpers work after installation.
* Added clearer runtime checks for `survdnn` when the LibTorch runtime is missing.
