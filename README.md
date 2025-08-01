# survalis



## survalis Learner Validation Checklist


* [X] `fit_aareg()`, `predict_aareg()`, `cv_survlearner()`, not tunable
* [X] `fit_aftgee()`, `predict_aftgee()`, `cv_survlearner()`, not tunable
* [X] `fit_bart()`, `predict_bart()`, `tune_bart()`
* [X] `fit_bnnsurv()`, `predict_bnnsurv()`, `tune_bnnsurv()`
* [X] `fit_cforest()`, `predict_cforest()`, `tune_cforest()`
* [X] `fit_coxaalen()`, `predict_coxaalen()`, `tune_coxaalen()`
* [X] `fit_cox()`, `predict_cox()`, `tune_coxph()`
* [X] `fit_flexsurv()`, `predict_flexsurv()`, `tune_flexsurv()`
* [X] `fit_glmnet()`, `predict_glmnet()`, `tune_glmnet()`
* [] `fit_mboost()`, `predict_mboost()`, `tune_mboost()`
* [] `fit_orsf()`, `predict_orsf()`, `tune_orsf()`
* [] `fit_ranger()`, `predict_ranger()`, `tune_ranger()`
* [] `fit_rfsrc()`, `predict_rfsrc()`, `tune_rfsrc()`
* [] `fit_rpart()`, `predict_rpart()`, `tune_rpart()`
* [] `fit_rstpm2()`, `predict_rstpm2()`, `tune_rstpm2()`
* [] `fit_selectcox()`, `predict_selectcox()`, `tune_selectcox()`
* [] `fit_survdnn()`, `predict_survdnn()`, `tune_survdnn()`
* [] `fit_survivalsvm()`, `predict_survivalsvm()`, `tune_survivalsvm()`
* [] `fit_survreg()`, `predict_survreg()`, `tune_survreg()`
* [] `fit_xgboost()`, `predict_xgboost()`, `tune_xgboost()`

---

### Remaining Tasks

* [ ] Validate `predict_*()` output: dimensions (`n x T`), column names (`t=...`), no NAs
* [ ] Validate `summary()` method for all learners
* [ ] Ensure all `tune_*()` functions use internal `cv_survlearner()` and return best model
* [ ] Verify all learners are registered in:
  * `cv_survlearner()`
  * `evaluate_survlearner()`
  * `benchmark_survlearners()`
  * Interpretation tools (`plot_shap`, `plot_pdp`, etc.)
* [ ] Add learners table to documentation: learner, engine, type, tuning
* [ ] Auto-test all learners via GitHub Actions



## NOTES 


- `aareg`

- fix max.time = max(data$time) and remove tuning. Otherwise this will gives inconsistent prediction horizons -> harder to compare within model




- `survivalsvm`






## `{survalis}` learners summary table 




| Learner        | Fit Function     | Predict Function     | Tuning Function   | Engine            | Type                       |
| -------------- | ---------------- | -------------------- | ----------------- | ----------------- | -------------------------- |
| coxph          | `fit_cox()`      | `predict_cox()`      | `tune_cox()`      | `survival`        | Semi-parametric (PH)       |
| aalen additive | `fit_aareg()`    | `predict_aareg()`    | `tune_aareg()`    | `survival`        | Additive hazards           |
| rsf            | `fit_rfsrc()`    | `predict_rfsrc()`    | `tune_rfsrc()`    | `randomForestSRC` | Tree ensemble              |
| cforest        | `fit_cforest()`  | `predict_cforest()`  | `tune_cforest()`  | `party`           | Tree ensemble              |
| rpart          | `fit_rpart()`    | `predict_rpart()`    | `tune_rpart()`    | `rpart`           | Single tree                |
| xgboost        | `fit_xgboost()`  | `predict_xgboost()`  | `tune_xgboost()`  | `xgboost`         | Boosted trees (AFT/Cox)    |
| gbm            | `fit_gbm()`      | `predict_gbm()`      | `tune_gbm()`      | `gbm`             | Boosted trees (Cox)        |
| flexsurv       | `fit_flexsurv()` | `predict_flexsurv()` | `tune_flexsurv()` | `flexsurv`        | Parametric (Weibull, etc.) |
| survreg        | `fit_survreg()`  | `predict_survreg()`  | `tune_survreg()`  | `survival/rms`    | Parametric AFT             |
| rstmp2         | `fit_rstpm2()`   | `predict_rstpm2()`   | `tune_rstpm2()`   | `rstpm2`          | Royston-Parmar (spline)    |
| mboost         | `fit_mboost()`   | `predict_mboost()`   | `tune_mboost()`   | `mboost`          | Boosted base-learners      |
| bart           | `fit_bart()`     | `predict_bart()`     | `tune_bart()`     | `BART`            | Bayesian Trees             |
| cox aalen      | `fit_coxaalen()` | `predict_coxaalen()` | `tune_coxaalen()` | `timereg`         | Cox + Additive             |
| survdnn        | `fit_survdnn()`  | `predict_survdnn()`  | `tune_survdnn()`  | `torch`           | Deep Neural Net (modular)  |





## Principales 

- No classes or S3 method dispatch

- Functions are pure: no side effects, no global state

- Output is plain data structures (lists + attributes)

- Prediction functions are composable



### CV design rational

- `feat: integrate evaluate_survlearner() with internal Surv construction`
- `feat: add unified cv_survlearner() for all engines`
- `refactor: drop evaluate_survmodel(), merge logic into CV loop`
- `feat: add support for multiple metrics in benchmark_survlearners()`
- `feat: add engine attribute to all fitted model objects`
- `feat: implement generic predict_survml() dispatcher`
- `refactor: standardize model object structure for mlsurv_model class`
- `feat: auto-handle status == event encoding in cv_survlearner()`



`cv_survlearner`: is concise CV engine for survival learners, designed to support all learners with `fit_*` and `predict_*.

Key design choices: 

- formula-first interface: aligns with base R modeling conventions and avoids manual specification of time/status
- flexible event handling: supports both raw status and filtered events
- automatic missing data handling: drops rows with NA in any variables used in the formula 
- tidy output: returns a consistent tiblle with metrics by fold, ready for summarization and plotting


### cv design metrics 

[] `evaluate_survlearner()`: provides a unified interface to compute survival probabilities at user-defined time points. 

[] fix warning in Brier related to empty folds in small sample size

[] build modular `cv_survlearner()`

[] build `benchmark_survleanrer()`  each learner tested via `purrr::pmap_dfr()`

[] build utilities for CV result summarization and best learner selection (`select_best()` + `cv_summary()` + `cv_plot()`)

[] handle edge case in Brier score when G(t) is NA or zero


### Learners implementation structure 



- `refactor: unify fit_* function signatures across learners`
- `feat: add consistent predict_* interface with times argument`
- `refactor: use Surv parsing logic to auto-detect time/status vars`
- `fix: ensure all predict_* return consistent t= colnames`

- `feat: add survlearners() registry for listing available learners`
- `feat: add select_best() for learner comparison by metric`
- `refactor: standardize colnames in prediction output (t= format)`
- `docs: document design rationale for standardized survival learners`
- `docs: add usage examples for fit_*, predict_*, and cv_* wrappers`



[] adds the `engine` attribute (`engine = <learner>`)
[] use consistent object structure (class = `"mlsurv_model"`)
[] keeps all naming conventions and time/status extraction compatible with the framework of this package 

### Rational behind learners 

- `aftgee`: semi-parametric AFT model via GEE

`fit_aftgee`
This learner implements a semiparametric accelerated failure time (AFT) model using generalized estimating equations via the aftgee package.
It allows robust estimation of survival models without specifying the full error distribution. This make it well suited for time-to-event with complex
structures, violation of PH assumption, and longitudinal or clustered survival data (although this implementation disables clustering for compatibility by fixing id = NULL)

`predict_aftgee`
The model output (TO CHECK) log-transformed survival time predictions, which are back-transformed to generate approximate survival probability
curves under a log-normal assumption. 


---

- `bart`

[] update the implementation (actual fit is to verbose)

---


- `mboost`

benchmarking through CV , `mboost::blackboost` may produce survival curves with time grid (`survFit$time`) that doesn't fully cover the user-specified `times`.
This causes inconsistencies and errors when evaluating metrics such as the c-index at fixed prediction times. 

To ensure consistency across folds and learners, we update `predict_mboost()` to: 

[] extend  the survival curve with the last know survival probability (flat extrapolation) if any `times` go beyond the model's range. This will avoid the metric computation failures due to missing time columns. 
[] add logic to handle times outside survFit$time via padding


---

- `cforest`

[] add conditional survival forest (`party:cforest`)
[] use interpolation to generate full survival matrices

--- 

- `ranger`

[] add random survival forest from `ranger`
[] ensure time alignment by nearest match 

---

- `orsf`

[] add oblique RF from aorsf 
[] use native `pred_horizon` and `pred_type = "surv"`


---

- `glmnet`

[] fit penalized cox (lasso/ridge) via via glmnet

[] manual estimated baseline hazard using `basehaz()`


---

- `flexsurv`

The `flexsurv` learner fits fully parametric survival models (e.g. weibull, log-normal, gompertz) using the `flexsurvreg()` function the `flexsurv` packagE. It provides smooth and flexible survival curves with closed-form expressions for prediction.

---



- `bnnsurv` 

Bagged Nearest Neighbors for survival analysis implements an ensemble of nearest-neighbor survival estimators using the `bnnSurvival` package. Each base learner uses a randomly seletect subset of features and training data to estimate individual survival curves.

This method is particularly well-suited for non-linear patterns and does not rely on PH assumptions. Predictions are interpolated over requested time points for full compatibility with the unified benchmarking pipeline.

## standardization of git commits fo


---

- `coxaalen` (additive-multiplicative hazards model)

This learner fits a cox-aalen model using `timereg::cox.aalen()` function, which combines both proportional hazards and additive effects. It allows for flexible modleing of time-varying and time-invariante covariates. 

_Preprocessing nots_: The initial `cox.aalen()` function requires explicit handling of covariates via `prop()` or `const()` to define how each covariates enters the model (as roportional hazards or additive effects). To avoid cruptic error like `Error: invalid type (NULL) for variable 'Z'`, our `fit_coxaalen()` wraps automatically all covariates with `prop()`, in order to ensure robust and concistent model fitting, enven if the user does not specify it manually. This ensure compatibility with bith simple and complex formulas.


--- 


- `rstpm2` (flexible parametric models)

This learner uses natural cubis splines on the log cumulative hazard scale for flexible survival modleing. 
It allows smoth survival function estimation and interpretable effects. The prediction method expands `newdata` over evaluation times and reshapes results into standard survival probability matrics.

_Note_: survival times ust be passed via the `newtime` argument of `predict()`, and expanded before calling the prediction method.


--- 


- `rpart` (survival tree)

This learner implement survival tree modeling using the `rpart` package with exponential survival (`method = "exp"`). 
For prediction, we approximate individual survival probabilities by assuming an exponential distribution for the predicted survival time per leaf node, then converting predicted node-level times into survival probabilities at user-specified evaluation times. 
This approach ensures uniform prediction output across all learners and maintains computational speed. Alternative KM-based similar to the approach implemented in `pec::pecRpart` may be considered in future updates if performance benchmarking or calibration tests suggest a meaningful benefit. 

_Note_: for now, we retain the exponentia-time approximation for simplicity and maintain consistency with other learners. 


---

- `rfsrc` (random survival forest)

Implements breiman's random forest methodology adapted for right-censored survival data.
Each tree is grown using a log-rank splitting rule, and survival probabilities are obtained by averaging ensemble kaplan-meier estimates across trees.
Predictions return interpolated survival probabilities at specified evaluation times. 


--- 

- `survivalsvm` (support vector machine for survival analysis)

This learner uses the `survivalsvm` package to fit a support vector machine model for survival data, using ranking or regression approach to estimate survival times. Since the model does not natively output survival probabilities, a workaround is applied by assuming a parametric form (exponential or weibull) to convert predicted survival times into survival probability curves. 

_Note_: this introduces an approximation that depends on the distributional assumption. The exponential model is the default, but a weibul shape can also be specified. 



---

- `xgboost` (extreme gradient boosting for survival)

This learner integrates the popular `xgboost` library to model survival data using either the Accelerated Failure Time (AFT) objective or the Cox PH hazards objective. It supports both right-censored outcomes and flexible distributions (`normal`, `logistic`, `extreme`) for AFT. 

Predictions are returned as survival probabilities at specified time points, computed by analytically transforming the model output margins usong the assumed distribution.

_Note_: The AFT model estimates the log-time to event, so survival probabilities are derived using parametric approximation (`1 - CDF(log(time))`). 



---

- `survreg` (parametric accelerated failure time AFT) 

This learner use the `rms::psm()` function to fit parametric accelerated failure time (AFT) survival models with flexible distributions (weibul, exponential, log-normal). 
Predictions returns survival probabilities at user specified time points. 




---


## Intepretability Methods

This package contains modular implementations of interpretablity techniques for survival models.
Each methos is organized in a separate `.R` script and includes: 

- `compute_*')`: function to generate interpretability data (e.g., partial dependence, ICE curves).
- `plot_*()`: functions to visualize the outputs using `ggplot2`.

### Currently implemented

[] pdp + ice
[] ale 
[] local surrogate model
[] tree surrogate
[] varimp
[] shap 
[] interaction 

---

### PDP: Partial Dependence Plots for survival models (methodological note)


- `feat: add compute_pdp() for PDP and ICE computation with survival predictions`
- `feat: add plot_pdp() for time-dependent and integrated PDP/ICE visualization`
- `feat: support both categorical and continuous features in compute_pdp()`
- `feat: implement integrated survival PDP using trapezoidal rule`
- `refactor: isolate compute_pdp() and plot_pdp() into separate interpretability files`
- `fix: handle single-row survival predictions in compute_pdp() robustly`
- `docs: add README note on survival-adapted PDP methodology`
- `test: add example usage of compute_pdp() with ranger on veteran dataset`
- `refactor: drop data.table dependency in compute_pdp() using base R reshape`


This implementation adapt partial dependence (pdp) and individual conditional expextation (ice) techniques to time-to-event data.

**Survival specific adaptation**

- The target fucntion is the survival probability $S(t|x)$ at a set of user defined time poits $t$.
- For each value `v` in the grid of a feature $X_j$, the model is used to compute $S(t | X[-j], X_j =v$ for all individuals. 

**PDP estimates**

- Mean survival probab accross individuals at each $(t, X_j = v)$.
- Displayed as time specific curves (facet) or summarized via numerical integratio over time ("integrated pdp").

**ICE estimates**

- individual-level survival predictions at each $(t, X_j = v)$.
- Useful to inspect heterogeneity in the effect of $X_j$ on survival.

**Integration strategy**

- For each feature value, we compute the integral of the predicted survival curve $S(t | X_j = v)$ across time. 
- This provides a summary score analogous to RMST, indicating the average survival time conditional on $X_j = v$. 
- Numerical integration is done using the trapezoidal rule. 

---

### Local surrogate interpretability (methodological note) 

- `feat: add penalized argument to compute_surrogate() with Lasso as default`
- `feat: support unpenalized local linear model when penalized = FALSE`
- `fix: robust exclusion of time/status using model formula`
- `feat: add exclude argument to manually specify variables to omit from surrogate fit`

This implementation provides a model-agnositic approach to locally explain survival predictions using a penalized linear surrogate model, inspired by LIME. 

For a given individual and target prediction time $t$, the methods: 

1) selects a local neighborhood from the training data using a gowerèbased distance kernel. 
2) predicts survival probabilities at time $t$ for each neighbor. 
3) encodes features relative to the invidual of interest (binary match for categorical, ra value for numeric).
4) Fits a sparse linear model (lasso or ols) to approximate the local behavior of the black-box model. 
5) Computes local feature effects as the product of the coefficients and the individual's encoded values. 

The output is a ranked list of the most influential features for the prediction, which enables interpretability inspection at the individual level. 

_Note_: options include control over the number of features selected $k$, distance metric, kernel, sharpness, and penalization.


### Tree surrogate interpretability (methodological note)

- `feat: add compute_tree_surrogate() to fit surrogate trees for survival models`

- `feat: add plot_tree_surrogate() for visualizing tree structure or split-count importance`

- `fix: ensure predict_ranger() handles one-time interpolation with consistent dimensions`

- `feat: support dynamic multi-time interpretation with tree surrogate summaries`

- `fix: robust check of prediction output in compute_tree_surrogate()`

- `feat: compute R² and split frequency as interpretability metrics`

This implementation supports time-specific interpretability using surrogate decision trees.

- `compute_tree_surrogate()`: Learns surrogate tree(s) that approximate survival predictions at specific time points.
- `plot_tree_surrogate()`: Visualizes decision trees (`type = "tree"`) or displays global feature importance via split counts (`type = "importance"`).

**Some main features**
- Compatible with any survival learner that returns a `data.frame` of survival probabilities.
- Supports both single and multiple time points (`times`) for dynamic model introspection.
- Computes R2 for surrogate fidelity and split frequency for feature importance.
- Relies on `rpart` (fitting) and `partykit` (visualization).


This approach builds a simple regression tree that approximates the predicted survival probability at a given time point $t$. The tree is fit using predicted survival probabilities as the response, and original covariates as predictors.

Let:
- $\hat{S}(t \mid X)$ be the predicted survival probability from a complex model at time $t$
- $X \in \mathbb{R}^p$ be the vector of input features

We define a surrogate model:

$$
\hat{S}(t \mid X) \approx f_{\text{tree}}(X)
$$

where $f_{\text{tree}}$ is a decision tree trained to predict $\hat{S}(t \mid X)$.

**Some advantages**

- Captures non-linear feature effects and interactions via tree splits
- Allows computation of surrogate $R^2$  to assess fidelity
- Global importance scores (via split frequency) can be aggregated across time

---

### Counterfactual explanation for survival models 

The compute_counterfactual() implements a method that identifies feature changes that would increase a patient's predicted survival probability at a specified time point. It supports both numeric and categorical features and uses a penalized gain objective: 

- survival gain: improvement in predicted survival probability from modifying a feature.
- change cost: size or effort of modifying the feature (absolute distance or 1 for categorical switch).
- penalized gain = survival gain - (penalty x change cost)

The function searches over a grid of possible feature values and selects the change that maximizes penalized gain. This enables actionable and interpretable "what-if" reasoning in clinical settings.

---

### Feature interaction diagnostics 

The `compute_interactions()` implements compute feature interaction strength using Freidman's H-statistic, extented to survival models with predicted survival probabilities at specific time points. 

Intraction stength is estimated using three approaches: 

- `1way`: each feature's total interaction with the rest (univariate H-statistic)
- `heatmap`: pairwise interactions visualized as a symmetric heatmap. 
- `time`: time-varying interaction strength per feature across time points. 

**Interpretation**

- A higher H-statistic (closer to 1) indicates stronger interaction (i.e., the feature's effect on survival depends more on other variables).
- Time-varying interaction results reveal how dependency structures evolve across observation period.


**Use cases**

- One can identify features with strong interactions to model them usong non-additive effects or interactions terms.
- Support time-varying effect modleing when interactions changes across time. 
- This may, in some extend, informs feature selection. 


---

## Tuning notes 


### Pattern for evaluating final model 

This pattern generalizes to all learners:

* `best_row <- res_* |> arrange(...) |> slice(1)`
* `final_model <- fit_*()` using those params
* Use `predict_*()` or `evaluate_survlearner()` as needed

**1. Extract the best row**

Assume you want to optimize for highest **C-index** (or lowest **IBS**, etc.):

```{r}
best_row <- res_cforest |>
  arrange(desc(cindex)) |>  # or arrange(ibs) if optimizing for IBS
  slice(1)
```


**2. Fit the final `cforest` model**

You pass the optimal parameters back into `fit_cforest()`:

```{r}
final_model_cforest <- fit_cforest(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  ntree = best_row$ntree,
  mtry = best_row$mtry,
  mincriterion = best_row$mincriterion,
  fraction = best_row$fraction
)
```

**3. Evaluate final model**

You can validate on the full dataset or on an external test set:

```{r}
evaluate_survlearner(
  model = final_model_cforest,
  metrics = c("cindex", "ibs", "iae", "ise"),
  times = c(100, 200, 300)
)
```
---

### No tuning for Aalen's Additive Model

We chose not to tune the max.time parameter for the Aalen additive hazards model (timereg::aalen) for methodological consistency. Unlike other learners, varying max.time changes the prediction horizon, leading to survival curves of different lengths and making metric comparisons unfair across models.

Additionally, timereg::aalen requires internal resampling for prediction, and tuning max.time can cause dimension mismatches in cross-validation outputs. To ensure stable and comparable evaluations, we fix max.time to the maximum follow-up time in the data and exclude it from hyperparameter tuning.


--

### No tuning for cox.aalen Model
We do not tune cox.aalen due to instability during cross-validation. 
Varying max.time alters the prediction grid, and low n.sim values lead to degenerate or missing predictions. As these parameters do not consistently improve predictive performance, we fix them to defaults and exclude this model from tuning.

---

### No tuning for survivalsvm 

We did not perform hyperparameter tuning for survivalsvm due to frequent optimization failures and unstable predictions across parameter combinations. The model was retained with default settings (gamma.mu = 0.1, kernel = "add_kernel", opt.meth = "ipop"), as tuning did not yield reliable improvements.

---

### No paramteric AFT model

We implemented the accelerated failure time (AFT) model using the rms::psm() function with a Weibull distribution. This learner provides stable and interpretable survival estimates, and supports direct estimation of survival probabilities at user-specified time points via survest(). No hyperparameter tuning was necessary as parametric models are fully specified once the distribution is chosen.

--------------------------------------------------------------------------------

# TODO


## STRUCTURE & METADATA

* [ ] Complete the `DESCRIPTION` file (title, description, version, authors, imports).
* [ ] Add `.onLoad()` in `zzz.R` to define package defaults.
* [ ] Ensure folder structure includes `R/`, `man/`, `tests/`, and `vignettes/`.
* [ ] Configure `.Rbuildignore` for non-package files (e.g., ZIPs, README\_files/).

---

## CORE IMPLEMENTATION

### Core infrastructure

* [ ] Review `fast_metrics.R` for speed and vectorization.
* [ ] Standardize input/output for `evaluate_survlearner()`.
* [ ] Modularize metric integration in `evaluate_survlearner()`.
* [ ] Add seed support to `cv_learner.R`.
* [ ] Add split validation and error checks to `cv_learner.R`.

---

### Learners

* [ ] Refactor `fit_all.R` to clean and centralize learner calls.
* [ ] Review and finalize `fit_coxaalen.R`.
* [ ] Review and finalize `fit_flexsurv.R`.
* [ ] Review and finalize `fit_mboost.R`.
* [ ] Review and finalize `fit_bart.R`.
* [ ] Review and finalize `fit_bnnsurv.R`.
* [ ] Review and finalize `fit_aareg.R`.
* [ ] Review and finalize `fit_aftgee.R`.
* [ ] Review and finalize `fit_cforest.R`.
* [ ] Review and finalize `fit_glmnet.R`.
* [ ] Ensure all `fit_*()` functions return `n x t` survival probability matrices.
* [ ] Add missing `predict_*()` functions for any `fit_*()` functions.
* [ ] Ensure all `predict_*()` functions handle both `data.frame` and `tibble`.

---

### Meta-learning

* [ ] Finalize `fit_metalearner()` function.
* [ ] Finalize `predict_metalearner()` function.
* [ ] Implement stacking and blending modes in `fit_metalearner()`.
* [ ] Add support for NNLS weights in meta-learner.
* [ ] Add support for ridge regression in meta-learner.
* [ ] Validate meta-learner predictions against base learners.

---

### Interpretability

* [ ] Finalize `compute_shap.R` implementation.
* [ ] Finalize `compute_ale.R` implementation.
* [ ] Finalize `compute_pdp.R` implementation.
* [ ] Finalize `compute_ceteris_paribus.R` implementation.
* [ ] Finalize `compute_varimp.R` implementation.
* [ ] Finalize `compute_interactions.R` implementation.
* [ ] Finalize `compute_surrogate.R` implementation.
* [ ] Finalize `compute_tree_surrogate.R` implementation.
* [ ] Implement `plot_shap()` for SHAP values.
* [ ] Implement `plot_ale()` for ALE curves.
* [ ] Implement `plot_pdp()` for PDPs.
* [ ] Implement `plot_varimp()` for feature importance.
* [ ] Implement `plot_ceteris_paribus()` for individual sensitivity.

---

### Causal Survival

* [ ] Finalize `compute_rmst.R` for RMST estimation.
* [ ] Finalize `compute_rmstdiff_tmle.R` for RMST diff with TMLE.
* [ ] Finalize `compute_rmstdiff_iptw.R` for RMST diff with IPTW.
* [ ] Finalize `compute_survdiff_gformula.R` for g-formula survival diff.
* [ ] Finalize `compute_survdiff_iptw.R` for IPTW survival diff.
* [ ] Finalize `compute_counterfactual.R` for counterfactual survival curves.
* [ ] Finalize `compute_iptw_weights.R` for stabilized weights.
* [ ] Finalize `compute_cate_rmst.R` for CATE via RMST.
* [ ] Finalize `compute_shap_cate.R` for CATE decomposition with SHAP.
* [ ] Implement `plot_shap_cate()` for causal SHAP visualization.

---

## DOCUMENTATION

* [ ] Write full Roxygen2 documentation for each exported function.
* [ ] Include `@param`, `@return`, `@examples`, and `@export` tags in each file.
* [ ] Tag internal helpers with `@keywords internal`.
* [ ] Create a package-level file `R/survalis-package.R`.
* [ ] Run `devtools::document()` after each module is finalized.

---

## TESTING

* [ ] Set up testing framework with `usethis::use_testthat()`.
* [ ] Write tests for all `fit_*()` functions.
* [ ] Write tests for all `predict_*()` functions.
* [ ] Write tests for `evaluate_survlearner()` and `cv_learner()`.
* [ ] Write tests for `fit_metalearner()` and `predict_metalearner()`.
* [ ] Write tests for interpretability functions (`compute_shap()`, etc.).
* [ ] Write tests for causal survival methods.
* [ ] Check test coverage with `covr::report()`.

---

## VIGNETTES & EXAMPLES

* [ ] Create vignette: “Getting Started with `survalis`”.
* [ ] Create vignette: “Model Evaluation & Meta-learners”.
* [ ] Create vignette: “Causal Survival & Interpretability”.
* [ ] Add examples to `README.Rmd` using minimal working code.

---

## 🧹 FINAL CHECKS & RELEASE

* [ ] Run `devtools::check()` and fix all warnings, errors, and notes.
* [ ] Run `R CMD check --as-cran` and verify.
* [ ] Run `goodpractice::gp()` for quality suggestions.
* [ ] Run `devtools::spell_check()` to fix typos.
* [ ] Finalize `NEWS.md` with full changelog for v0.1.0.
* [ ] Tag GitHub release as `v0.1.0`.

---

## VERSIONING & PUBLICATION

* [ ] Start development version as `0.0.0.9000`.
* [ ] Mark first stable CRAN release as `0.1.0`.
* [ ] Add `Lifecycle: experimental` badge to README.
* [ ] Tag `v1.0.0` once stable and CRAN-released.
* [ ] Submit JOSS paper after CRAN release is live.

---

## PROOFREADING, CONSISTENCY & FINAL CHECKS

* [ ] Review all function names and ensure consistent naming conventions.
* [ ] Check for duplicated logic or redundant functions across files.
* [ ] Ensure consistent argument names across all `fit_*()` and `predict_*()` functions.
* [ ] Verify that all `predict_*()` functions return properly formatted survival matrices.
* [ ] Cross-check that all `compute_*()` functions work with standardized inputs: `model`, `newdata`, `times`.
* [ ] Verify that all `plot_*()` functions accept output from the corresponding `compute_*()` functions.
* [ ] Manually run all example code in Roxygen `@examples` to ensure they execute cleanly.
* [ ] Rebuild documentation with `devtools::document()` and inspect generated `.Rd` files.
* [ ] Run a full test suite with `devtools::test()` and ensure no failures.
* [ ] Run a final `devtools::check()` and resolve all errors, warnings, and notes.
* [ ] Verify vignette builds with `devtools::build_vignettes()` and `pkgdown::build_site()`.
* [ ] Confirm CRAN readiness with `R CMD check --as-cran` on local machine.
* [ ] Perform a final spell check using `devtools::spell_check()`.
* [ ] Review README, vignettes, and NEWS.md for accuracy, formatting, and coherence.
* [ ] Cross-check all internal functions are unexported and tagged `@keywords internal`.
