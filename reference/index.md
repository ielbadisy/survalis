# Package index

## Workflow verbs

The high-level interface: fit a learner, evaluate it, compare learners,
tune, interpret and estimate contrasts from one formula and one dataset.

- [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) : Fit a
  survival learner
- [`predict(`*`<survalis_fit>`*`)`](https://elbadisyimad.com/survalis/reference/predict.survalis_fit.md)
  : Predict from a fitted survival learner
- [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
  : Resampled performance of a survival learner
- [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  : Compare several survival learners on shared resampling folds
- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) : Tune
  a survival learner over a hyperparameter grid
- [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
  : Interpret a fitted survival learner
- [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)
  : Standardised causal-style contrast for a binary treatment
- [`as.data.frame(`*`<survalis_eval>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_tune>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_compare>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  [`as.data.frame(`*`<survalis_estimate>`*`)`](https://elbadisyimad.com/survalis/reference/survalis-tidiers.md)
  : Tidy data frames from verb result objects
- [`survalis`](https://elbadisyimad.com/survalis/reference/survalis-package.md)
  [`survalis-package`](https://elbadisyimad.com/survalis/reference/survalis-package.md)
  : survalis: Interpretable Survival Machine Learning

## Resampling

Split schemes used by
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
and [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md).

- [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  : Resampling specifications for evaluate(), tune() and compare()

## Learners

Low-level fit and predict pairs, one per engine. Every `predict_*()`
returns a survival matrix (rows are subjects, columns are times).

- [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md)
  : List survival learners available in survalis

- [`summary(`*`<mlsurv_model>`*`)`](https://elbadisyimad.com/survalis/reference/summary.mlsurv_model.md)
  :

  Summarize an `mlsurv_model`

## Tuning

Grid search per learner.

- [`list_tunable_survlearners()`](https://elbadisyimad.com/survalis/reference/list_tunable_survlearners.md)
  : List tunable survival learners

## Cross-validation and benchmarking

- [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
  :

  Cross-Validate a Survival Learner (fold-mapped with `fmapn`)

- [`cv_survmetalearner()`](https://elbadisyimad.com/survalis/reference/cv_survmetalearner.md)
  : Cross‑Validate a Stacked Survival Meta‑Learner

- [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md)
  : Summarize Cross-Validation Results

- [`cv_plot()`](https://elbadisyimad.com/survalis/reference/cv_plot.md)
  : Boxplot of Cross-Validation Metric Distributions

- [`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md)
  : Benchmark Tuned Survival Learners with Nested Cross-Validation

- [`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md)
  : Summarise Benchmark Results (Mean SD with Wald CI)

- [`summarize_benchmark_results()`](https://elbadisyimad.com/survalis/reference/summarize_benchmark_results.md)
  : Compact Table of Mean SD by Learner and Metric

- [`best_survlearner()`](https://elbadisyimad.com/survalis/reference/best_survlearner.md)
  : Select the Best Survival Learner by a Given Metric

- [`plot_benchmark()`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md)
  : Plot Benchmark Distributions Across Learners

- [`plot_survmetalearner_weights()`](https://elbadisyimad.com/survalis/reference/plot_survmetalearner_weights.md)
  : Plot Time‑Varying Stacking Weights

## Performance metrics

Discrimination, calibration and overall accuracy on a survival matrix.

- [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md)
  : Score a Fitted Survival Model on Its Training Data
- [`list_metrics()`](https://elbadisyimad.com/survalis/reference/list_metrics.md)
  : List Available Evaluation Metrics
- [`cindex_survmat()`](https://elbadisyimad.com/survalis/reference/cindex_survmat.md)
  : Concordance Index from a Survival-Probability Matrix
- [`brier()`](https://elbadisyimad.com/survalis/reference/brier.md) :
  Brier Score with IPCW for a Single Time Point
- [`ibs_survmat()`](https://elbadisyimad.com/survalis/reference/ibs_survmat.md)
  : Integrated Brier Score (Discrete Integration)
- [`iae_survmat()`](https://elbadisyimad.com/survalis/reference/iae_survmat.md)
  : Integrated Absolute Error Against Kaplan-Meier
- [`ise_survmat()`](https://elbadisyimad.com/survalis/reference/ise_survmat.md)
  : Integrated Squared Error Against Kaplan-Meier
- [`auc_survmat()`](https://elbadisyimad.com/survalis/reference/auc_survmat.md)
  : Time-Dependent AUC from a Survival-Probability Matrix
- [`roc_survmat()`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  [`plot(`*`<roc_survmat>`*`)`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  [`plot_roc()`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  : Time-Dependent ROC Curve from a Survival-Probability Matrix
- [`timeroc_survmat()`](https://elbadisyimad.com/survalis/reference/timeroc_survmat.md)
  : Time-Dependent ROC/AUC Curve from a Survival-Probability Matrix
- [`ece_survmat()`](https://elbadisyimad.com/survalis/reference/ece_survmat.md)
  : Expected Calibration Error (ECE) at a Single Time Point
- [`dca_survmat()`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  [`plot(`*`<dca_survmat>`*`)`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  [`plot_dca()`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  : Decision Curve Analysis for Right-Censored Survival Outcomes
- [`compute_calibration()`](https://elbadisyimad.com/survalis/reference/compute_calibration.md)
  : Calibration of Survival Predictions at a Single Time Point

## Interpretability

Model-agnostic explanation methods and their plots.

- [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md)
  : List interpretability methods available in survalis
- [`compute_ale()`](https://elbadisyimad.com/survalis/reference/compute_ale.md)
  : Accumulated Local Effects (ALE) for Survival Models
- [`compute_calibration()`](https://elbadisyimad.com/survalis/reference/compute_calibration.md)
  : Calibration of Survival Predictions at a Single Time Point
- [`compute_counterfactual()`](https://elbadisyimad.com/survalis/reference/compute_counterfactual.md)
  : Compute individual counterfactual changes to increase survival
- [`compute_interactions()`](https://elbadisyimad.com/survalis/reference/compute_interactions.md)
  : Compute Feature Interactions for Survival Predictions
- [`compute_pdp()`](https://elbadisyimad.com/survalis/reference/compute_pdp.md)
  : Partial Dependence and ICE for Survival Predictions
- [`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md)
  : Compute local SHAP-like contributions for survival predictions
- [`compute_shap_matrix()`](https://elbadisyimad.com/survalis/reference/compute_shap_matrix.md)
  : Time-Integrated SHAP Values for a Set of Observations
- [`compute_surrogate()`](https://elbadisyimad.com/survalis/reference/compute_surrogate.md)
  : Local Surrogate Explanation for Survival Predictions (LIME-style)
- [`compute_survcluster()`](https://elbadisyimad.com/survalis/reference/compute_survcluster.md)
  : Cluster Predicted Survival Curves
- [`compute_tree_surrogate()`](https://elbadisyimad.com/survalis/reference/compute_tree_surrogate.md)
  : Compute Tree-Based Surrogate Model for Survival Predictions
- [`compute_varimp()`](https://elbadisyimad.com/survalis/reference/compute_varimp.md)
  : Permutation variable importance for survival models
- [`plot_varimp()`](https://elbadisyimad.com/survalis/reference/plot_varimp.md)
  : Plot Permutation Variable Importance
- [`plot_pdp()`](https://elbadisyimad.com/survalis/reference/plot_pdp.md)
  : Plot PDP/ICE Curves for Survival Models
- [`plot_ale()`](https://elbadisyimad.com/survalis/reference/plot_ale.md)
  : Plot ALE Curves for Survival Models
- [`plot_shap()`](https://elbadisyimad.com/survalis/reference/plot_shap.md)
  : Plot SHAP-like contributions for survival models
- [`plot_shap_beeswarm()`](https://elbadisyimad.com/survalis/reference/plot_shap_beeswarm.md)
  : SHAP Beeswarm Plot (Time-Integrated)
- [`plot_interactions()`](https://elbadisyimad.com/survalis/reference/plot_interactions.md)
  : Plot Interaction Strengths for Survival Models
- [`plot_counterfactual()`](https://elbadisyimad.com/survalis/reference/plot_counterfactual.md)
  : Plot Counterfactual Recommendations
- [`plot_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_surrogate.md)
  : Plot Local Surrogate Explanation
- [`plot_tree_surrogate()`](https://elbadisyimad.com/survalis/reference/plot_tree_surrogate.md)
  : Plot Tree-Based Surrogate Models or Feature Importances
- [`plot_survcluster()`](https://elbadisyimad.com/survalis/reference/plot_survcluster.md)
  : Plot Clustered Survival Curves

## Survival matrices and estimands

Convert predicted survival into hazards, quantiles and RMST.

- [`new_survmat()`](https://elbadisyimad.com/survalis/reference/survmat.md)
  [`validate_survmat()`](https://elbadisyimad.com/survalis/reference/survmat.md)
  [`as_survmat()`](https://elbadisyimad.com/survalis/reference/survmat.md)
  : The survmat prediction contract
- [`survmat_times()`](https://elbadisyimad.com/survalis/reference/survmat_times.md)
  : Evaluation grid of a survmat
- [`survmat_to_chf()`](https://elbadisyimad.com/survalis/reference/survmat_to_chf.md)
  : Convert a survival-probability matrix (survmat) to cumulative hazard
- [`survmat_to_haz()`](https://elbadisyimad.com/survalis/reference/survmat_to_haz.md)
  : Convert a survival-probability matrix (survmat) to hazards on a time
  grid
- [`survmat_to_quantile()`](https://elbadisyimad.com/survalis/reference/survmat_to_quantile.md)
  : Compute a survival-time quantile from a survival-probability matrix
  (survmat) with grid-based approach
- [`survmat_to_rmst()`](https://elbadisyimad.com/survalis/reference/survmat_to_rmst.md)
  : Compute restricted mean survival time (RMST) from a
  survival-probability matrix (survmat)
- [`estimate_rmst()`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md)
  : Standardized RMST Estimate or Causal-Style RMST Contrast via
  G-Computation
- [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md)
  : Default evaluation-time grid
- [`screen_fdr()`](https://elbadisyimad.com/survalis/reference/screen_fdr.md)
  : FDR-Based Variable Screening via Univariate Cox Models

## Plotting

- [`plot_survcurve()`](https://elbadisyimad.com/survalis/reference/plot_survcurve.md)
  : Kaplan-Meier Survival Curve with Risk Table (survminer-Inspired)

- [`plot_survmat()`](https://elbadisyimad.com/survalis/reference/plot_survmat.md)
  : Plot Predicted Survival Curves from a survmat

- [`plot_calibration()`](https://elbadisyimad.com/survalis/reference/plot_calibration.md)
  : Plot Calibration Curve for Survival Predictions

- [`roc_survmat()`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  [`plot(`*`<roc_survmat>`*`)`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  [`plot_roc()`](https://elbadisyimad.com/survalis/reference/roc_survmat.md)
  : Time-Dependent ROC Curve from a Survival-Probability Matrix

- [`dca_survmat()`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  [`plot(`*`<dca_survmat>`*`)`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  [`plot_dca()`](https://elbadisyimad.com/survalis/reference/dca_survmat.md)
  : Decision Curve Analysis for Right-Censored Survival Outcomes

- [`plot_rmst()`](https://elbadisyimad.com/survalis/reference/plot_rmst.md)
  : RMST Estimate or Contrast Curve Over a Grid of Restriction Times

- [`theme_survalis()`](https://elbadisyimad.com/survalis/reference/theme_survalis.md)
  : A Consistent ggplot2 Theme for survalis Figures

- [`scale_color_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md)
  [`scale_fill_survalis()`](https://elbadisyimad.com/survalis/reference/survalis-scales.md)
  :

  Discrete Color/Fill Scales Matching `theme_survalis`

## Low-level interface

- [`survalis-lowlevel`](https://elbadisyimad.com/survalis/reference/survalis-lowlevel.md)
  : Low-level functions in survalis

## Data and re-exports

- [`veteran`](https://elbadisyimad.com/survalis/reference/veteran.md) :
  Veteran's Administration Lung Cancer Trial Data
- [`reexports`](https://elbadisyimad.com/survalis/reference/reexports.md)
  [`Surv`](https://elbadisyimad.com/survalis/reference/reexports.md) :
  Re-export Surv from survival
