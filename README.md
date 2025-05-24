# isurvml


## Principales 

- No classes or S3 method dispatch

- Functions are pure: no side effects, no global state

- Output is plain data structures (lists + attributes)

- Prediction functions are composable





### CV design rational


`cv_survlearner`: is concise CV engine for survival learners, designed to support all learners with `fit_*` and `predict_*.

Key design choices: 

- formula-first interface: aligns with base R modeling conventions and avoids manual specification of time/status
- flexible event handling: supports both raw status and filtered events
- automatic missing data handling: drops rows with NA in any variables used in the formula 
- tidy output: returns a consistent tiblle with metrics by fold, ready for summarization and plotting


### cv design metrics 

[X] `evaluate_survlearner()`: provides a unified interface to compute survival probabilities at user-defined time points. 

[] fix warning in Brier related to empty folds in small sample size

[] build modular `cv_survlearner()`

[] build `benchmark_survleanrer()`  each learner tested via `purrr::pmap_dfr()`

[] build utilities for CV result summarization and best learner selection (`select_best()` + `cv_summary()` + `cv_plot()`)

[] handle edge case in Brier score when G(t) is NA or zero


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

[X] extend  the survival curve with the last know survival probability (flat extrapolation) if any `times` go beyond the model's range. This will avoid the metric computation failures due to missing time columns. 
[X] add logic to handle times outside survFit$time via padding


---

- `cforest`

[X] add conditional survival forest (`party:cforest`)
[X] use interpolation to generate full survival matrices

--- 

- `ranger`

[X] add random survival forest from `ranger`
[X] ensure time alignment by nearest match 

---

- `orsf`

[X] add oblique RF from aorsf 
[X] use native `pred_horizon` and `pred_type = "surv"`


---

- `glmnet`

[X] fit penalized cox (lasso/ridge) via via glmnet

[X] manual estimated baseline hazard using `basehaz()`


---

- `flexsurv`

The `flexsurv` learner fits fully parametric survival models (e.g. weibull, log-normal, gompertz) using the `flexsurvreg()` function the `flexsurv` packagE. It provides smooth and flexible survival curves with closed-form expressions for prediction.

---

## standardization of git commits fo


**Model API**

- `refactor: unify fit_* function signatures across learners`
- `feat: add consistent predict_* interface with times argument`
- `refactor: use Surv parsing logic to auto-detect time/status vars`
- `fix: ensure all predict_* return consistent t= colnames`

**Learners**

- `feat: add engine attribute to all fitted model objects`
- `feat: implement generic predict_survml() dispatcher`
- `refactor: standardize model object structure for mlsurv_model class`
- `feat: auto-handle status == event encoding in cv_survlearner()`

- **Evaluation and CV**

- `feat: integrate evaluate_survlearner() with internal Surv construction`
- `feat: add unified cv_survlearner() for all engines`
- `refactor: drop evaluate_survmodel(), merge logic into CV loop`
- `feat: add support for multiple metrics in benchmark_survlearners()`


**Usability**

- `feat: add survlearners() registry for listing available learners`
- `feat: add select_best() for learner comparison by metric`
- `refactor: standardize colnames in prediction output (t= format)`
- `docs: document design rationale for standardized survival learners`
- `docs: add usage examples for fit_*, predict_*, and cv_* wrappers`





