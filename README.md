
## Implementation design and core philosophies

The `survalis` framework is built on a foundational principle of _standardization and modularity_, enabling seamless comparison, interpretation, and benchmarking of a wide range of survival learning algorithms. All survival learners are implemented using a consistent trio of functions: `fit_*()`, `predict_*()`, and `tune_*()`, with standardized inputs and outputs. Each fitted model is wrapped in a minimal `mlsurv_model` object containing the original formula, training data, engine name, and learner identifier. This approach ensures predictability, composability, and ease of downstream evaluation.

We deliberately avoid S3 class hierarchies or method dispatch. Instead, we use plain R functions that return standardized lists with required metadata. This choice enforces purity (no side effects), simplifies debugging, and enhances transparency. Prediction functions always return a plain `data.frame` with numeric survival probabilities, formatted with column names of the form `t=100`, `t=200`, etc., ensuring compatibility with all metric computation routines. 



## Learner-specific design choices

Some learners require special handling to ensure compatibility and stability. For instance, we fix `max.time = max(data$time)` in the Aalen additive hazards model (`aareg`) to prevent inconsistent prediction horizons across folds or resamples. Similarly, for `cox.aalen` models from `timereg`, we automatically wrap all covariates in `prop()` to ensure consistent proportional hazard modeling and prevent obscure errors related to missing `Z` variables.

The `survivalsvm` model, which lacks native survival probability outputs, is handled via _parametric approximation_: we convert predicted survival times into survival probabilities under an exponential or Weibull assumption. While approximate, this ensures output consistency for comparison across learners. Due to frequent instability, no hyperparameter tuning is performed for this learner.

For tree-based models like `mboost` or `rfsrc`, we implement custom interpolation and padding logic to ensure predicted survival probabilities align exactly with user-specified time grids. This is essential for fair metric computation and avoids runtime failures in validation pipelines (e.g., when evaluating c-index or IBS).


## Cross-Validation and evaluation framework

The central evaluation routine in `survalis` is `cv_survlearner()`, which performs cross-validation of survival learners in a fully standardized manner. Key design principles include a formula-first interface, automatic handling of missing data, flexible support for status/event coding, and tidy outputs suitable for downstream analysis and plotting. All model training, prediction, and evaluation are encapsulated within the CV loop, and parallelism is centralized in this function to avoid nested concurrency issues.

Parallel execution is limited to `cv_survlearner()` to prevent CPU oversubscription and to provide users with a clear mental model of where computation is concentrated. Tuning functions (`tune_*()`) remain sequential and rely on cross-validation results to guide hyperparameter selection.



## Prediction interface

The `predict_*()` functions follow a unified contract: given a trained `mlsurv_model`, a new dataset, and a vector of prediction times, they return a rectangular `data.frame` with predicted survival probabilities. The design ensures that row = individual, column = time point regardless of the model. Special handling is implemented for learners that do not support arbitrary time prediction (e.g., `mboost`, `rpart`, or `bnnsurv`) through interpolation or extrapolation strategies to match requested times.

This standardization is critical for enabling model-agnostic interpretability, benchmarking, and diagnostic tools. It ensures that interpretability methods like PDP, SHAP, and surrogate models operate on a consistent prediction interface.


## Interpretability method architecture

Each interpretability method in `survalis` is implemented via a `compute_*()` + `plot_*()` pairing, following a clear separation of concerns. Methods like partial dependence (`compute_pdp()`), accumulated local effects (`compute_ale()`), and feature importance (`compute_varimp()`) rely on the unified prediction interface. Advanced methods like local surrogate models and tree-based surrogates build explanations at either individual or global levels.

The design supports both time-specific and integrated views of interpretability, leveraging numerical integration (e.g., trapezoidal rule) to summarize survival curves into scalar metrics when needed (e.g., for RMST-like interpretations in PDPs).


## Benchmarking and model selection

For consistent learner comparison, we provide `benchmark_survlearners()`, which internally relies on `cv_survlearner()` and aggregates performance across a list of learners using `purrr::pmap_dfr()`. This enables comprehensive, reproducible benchmarking using multiple metrics (C-index, IBS, IAE, ISE). Supplementary tools like `select_best()`, `cv_summary()`, and `cv_plot()` help users inspect and summarize results.

Hyperparameter tuning follows a consistent pattern: extract the best row based on a metric, refit the final model using `fit_*()`, and evaluate using `evaluate_survlearner()` or `score_survmodel()`. Learners that are unstable or fully parametric (e.g., `aareg`, `cox.aalen`, `survreg`) are excluded from tuning, as tuning does not yield meaningful improvements and can introduce unnecessary complexity or inconsistencies.


## Remarks

Overall, the `survalis` design philosophy emphasizes clarity and modularity. Each component, from learners and evaluators to interpretability tools—follows a unified interface. This enables practitioners and researchers to focus on model logic and decision-making, not debugging infrastructure mismatches. With this system, adding a new learner or interpretability method becomes a matter of implementing a minimal set of well-defined functions, backed by robust and reusable tooling.
