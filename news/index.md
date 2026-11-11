# Changelog

## survalis 2.3.0

- Exported names are glued lowercase, like the verbs:
  [`listlearners()`](https://elbadisyimad.com/survalis/reference/listlearners.md),
  [`listmetrics()`](https://elbadisyimad.com/survalis/reference/listmetrics.md),
  [`listmethods()`](https://elbadisyimad.com/survalis/reference/listmethods.md)
  and
  [`groupcv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  replace
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  [`list_metrics()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md)
  and
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md).
  The old names still work with a deprecation warning and will be
  removed in the next release.

### Bug fixes

- `coxph` predicts with `strata()` in the formula (one curve per row).
- `aftgee` fits training data with missing covariates; its default
  prediction times come from the training data.
- `ranger` passes `case.weights` to the engine.
- `blackboost` keeps mboost’s unrestricted tree tests; trees were
  stopped early by a Bonferroni test.
- `survsvm` stores the correct outcome column names.
- `densemlp` reads its exact time grid; large times were rounded to four
  digits.
- `survmetalearner` uses equal weights when NNLS returns all zeros
  instead of producing NaN.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) of a tuning
  result draws one line per value of the second tuned parameter instead
  of joining all configurations.

## survalis 2.2.2

- `"aftgee"` predictions now use the semiparametric baseline of the AFT
  model: the Kaplan-Meier estimate of the censored residuals log(T) -
  x’beta, read at log(t) - x’beta. They previously assumed a log-normal
  error with scale fixed at 1, which miscalibrated survival: on the
  Rotterdam breast cancer data, mean predicted 5-year survival fell from
  0.47 to 0.44 against a Kaplan-Meier estimate of 0.44.
- `"bart"` predictions between the learner’s coarse time-grid points are
  now interpolated on the log-survival scale, as for `"densemlp"` in
  2.2.1, instead of holding the value from the start of the interval.

## survalis 2.2.1

- `"densemlp"` predictions between the learner’s coarse time-grid points
  are now interpolated on the log-survival scale (a constant hazard
  within each bin) instead of holding the value from the start of the
  bin. With the default 10 bins the old lookup overstated survival by up
  to one bin: on the Rotterdam breast cancer data, mean predicted 5-year
  survival fell from 0.54 to 0.42 against a Kaplan-Meier estimate of
  0.44.

## survalis 2.2.0

- Removed functions that the verbs replace: `benchmark()`,
  `benchmark_default_survlearners()`, `benchmark_tuned_survlearners()`,
  `summarise_benchmark()`, `summarize_benchmark_results()`,
  `best_survlearner()` and `plot_benchmark()` (use
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)),
  `estimate_rmst()` and `plot_rmst()` (use
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)),
  `cv_survmetalearner()` (use `evaluate(model = "survmetalearner")`),
  `plot_survmetalearner_weights()`, `compute_shap_mean()` (use
  `interpret(method = "shap")`), `plot_survcurve()` (use
  [`kmplot()`](https://elbadisyimad.com/survalis/reference/kmplot.md))
  and `screen_fdr()` (use the `"selectcox"` or `"glmnet"` learners).

- New
  [`kmplot()`](https://elbadisyimad.com/survalis/reference/kmplot.md):
  Kaplan-Meier curves in a publication style, with censoring marks, an
  optional log-rank p-value, optional confidence bands and a
  number-at-risk table coloured by group. `type` draws the survival
  function, the risk function (one minus survival), the cumulative
  hazard, or the log(-log) survival against log time for checking
  proportional hazards. A factor status gives the Aalen-Johansen
  cumulative incidence of the chosen `event` under competing risks. A
  fitted `model` adds its curves, dashed: the mean prediction within
  each group (`model.curve = "predicted"`, a calibration check) or
  covariate-adjusted curves by G-computation
  (`model.curve = "adjusted"`). Panels combine with ‘patchwork’
  operators. survalis now requires ‘ggplot2’ 3.5.0 or later.

### Bug fixes

These change results; rerun analyses made with earlier versions.

- Outcome handling:
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) arguments must
  now be column names. A computed time or status
  (`Surv(time / 365, status)`) was silently read from the raw column and
  is now an error. A 1/2 or logical status (as in
  [`survival::lung`](https://rdrr.io/pkg/survival/man/lung.html)) and a
  `status == k` recoding now reach every learner and metric as a 0/1
  event. [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md)
  drops incomplete rows once, and `~ .` works in
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)
  and the explanation methods.

- Metrics: the C-index now equals
  [`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html),
  and the AUC weights cases by the censoring survival just before their
  event time, as in ‘timeROC’. The IBS no longer depends on the order of
  `times`; IAE and ISE integrate over the evaluation window only. Metric
  values are stored at full precision (rounding happens only when
  printing), the standard error counts only non-missing splits, and
  unknown metric names, or single-time metrics given several times, are
  rejected before any model is fitted.

- Resampling: each split gets its own seed, so
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  and [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md)
  give the same values on repeated runs and whatever `ncores` is.
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) keeps
  engine arguments passed in `...` in the best spec and rejects an
  argument that is both tuned and fixed.

- Cox-type predictions: `"coxph"` now uses
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html),
  and the penalized Cox learners use a Breslow baseline on the training
  linear predictor. The earlier baseline was too large when the mean
  linear predictor was far from 0.

- Explanations: `interpret(method = "shap")` uses a permutation
  estimator whose values add up to the prediction gap; `"ale"` follows
  Apley and Zhu (2020), using only the observations inside each bin;
  `"interactions"` computes the centred Friedman-Popescu H statistic;
  `"varimp"` scores the fixed model on permuted data and reports a
  positive importance when permuting a feature makes predictions worse,
  whatever the metric’s direction. Calibration uses `newdata` when given
  and handles tied predictions. `"shap"`, `"counterfactual"` and
  `"surrogate"` require a one-row `newdata`, and the explanation methods
  use only the model’s own predictors.

- Learners:

  - Tree, forest and network learners (`"rsf"`, `"ranger"`, `"cforest"`,
    `"bart"`, `"blackboost"`, `"densemlp"`) read their curves as step
    functions, equal to 1 before the first time.
  - `"rpart"` uses the Kaplan-Meier curve of each leaf.
  - `"nbsurv"` predicts each row on its own, avoiding an engine bug in
    batch prediction.
  - `"rsf"` follows [`set.seed()`](https://rdrr.io/r/base/Random.html).
  - `"survmetalearner"` stacks on subjects with known status at each
    landmark, with censoring weights.
  - `"xgboost"` with `objective = "survival:cox"` codes censored times
    as xgboost expects and uses a Breslow baseline. With
    `"survival:aft"` it starts from the geometric mean time.
  - `"bart"` passes `...` to
    [`BART::surv.bart()`](https://rdrr.io/pkg/BART/man/surv.bart.html).
  - `"selectcox"` carries survival forward past the last follow-up.
  - The design-matrix learners (`"glmnet"`, `"coxboost"`, `"xgboost"`,
    `"aalen"`, `"aftgee"`) keep the training factor levels, so a single
    row or a test fold missing a level is coded correctly, and `newdata`
    no longer needs the outcome columns.

- `kmplot(model = )` accepts a model fitted on a rescaled time (days
  versus years) and warns when the model does not use the grouping
  variable.

## survalis 2.0.0

- survalis is now built around seven verbs:
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md),
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md),
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
  [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
  and
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md).
  The exported API is these verbs plus the resampling specs
  ([`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  [`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md)),
  the lookups
  ([`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  [`list_metrics()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md)),
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) and the
  `veteran` data. The granular layer is no longer exported: `fit_*()`,
  `predict_*()`, `tune_*()`, `cv_survlearner()`, `cv_summary()`,
  `cv_plot()`, `score_survmodel()`, `compute_*()`, `plot_*()`, the
  `*_survmat()` metrics, `brier()`, `survmat_to_*()`, `default_times()`,
  `list_tunable_survlearners()`, `theme_survalis()` and the colour
  scales.
  [`vignette("from-1.0")`](https://elbadisyimad.com/survalis/articles/from-1.0.md)
  maps each one to its verb.

- Table operations inside survalis (reshaping, grouped summaries, row
  binding, sorting and merging) now use ‘basetable’ instead of base R.
  Result tables stored in fitted objects
  ([`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
  [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md),
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md))
  are basetable objects;
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) on a
  survalis result still returns a plain `data.frame`.

- New decision-weighted scores
  ([`dwscore()`](https://elbadisyimad.com/survalis/reference/dwscore.md)):
  the loss of the rule “act if predicted risk exceeds c”, averaged over
  risk thresholds c with a stated weight (threshold range, SurvAcc,
  Brier, logarithmic or Beta), with IPCW for censoring. Each score is a
  weighted area under the net benefit curve. The metric ids `"survacc"`
  and `"range:<lower>-<upper>"` (for example `"range:0.05-0.20"`) work
  in
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  and [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md).

- The `"brier"` and `"ibs"` metrics now use the same censoring weights
  as
  [`dwscore()`](https://elbadisyimad.com/survalis/reference/dwscore.md)
  and the usual IPCW Brier score (Graf et al. 1999): events at or before
  `t` are weighted by `1 / G(T-)`. Values can change in the fourth
  decimal.

- The decision curve (`interpret(method = "dca")`) now uses the same
  censoring weights, treats a patient when predicted risk is above the
  threshold (not equal to it), and computes “treat all” from the same
  weights instead of the Kaplan-Meier risk. Each
  `"range:<lower>-<upper>"` loss then equals the weighted area between
  the treat-all and model curves exactly.

- Custom learners: `model` in every verb also accepts
  `list(fit = function(formula, data, ...), predict = function(object, newdata, times))`,
  where `predict` returns an `nrow(newdata) x length(times)` matrix of
  survival probabilities.
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  accepts a list mixing ids and custom learners.

- `evaluate(fit, newdata = test)` scores a fitted learner on new data
  (external or temporal validation), and
  `evaluate(S, y = Surv(time, status))` scores a survival matrix
  produced anywhere.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) also
  works on a
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) result
  (one row: `model`, `engine`, `formula`, `n`, `events`; it was an
  error), and `as.data.frame(predict(...), format = "long")` gives one
  row per observation and time (`id`, `time`, `value`). The
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)
  result and its summary name the interval bounds `conf_low` and
  `conf_high` (were `conf.low` and `conf.high`).

- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) works
  for any learner when you supply `grid`: configurations are scored on
  shared splits under any resampling spec. Learners listed by
  `list_survlearners(has_tune = TRUE)` keep their built-in default
  grids.

- A `grid` passed to
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) no
  longer has to name every parameter of a learner’s built-in tuner. A
  grid that left one out (for example
  `grid = data.frame(cp = c(0.001, 0.01))` for `"rpart"`) stopped with
  ‘argument “maxdepth” is missing’. Parameters not in the grid now take
  the learner’s defaults: such a grid, or one with extra engine
  arguments, is scored on shared splits like the grid of any other
  learner, and `...` supplies fixed arguments to every configuration. A
  grid with exactly the built-in columns still uses the built-in tuner.
  The `"survsvm"` and `"glmnet"` tuners also take a data frame grid row
  by row instead of crossing its columns.

- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) with a
  built-in default grid returned no configurations for the `"glmnet"`,
  `"bart"`, `"blackboost"`, `"bnnsurv"`, `"cforest"`, `"flexsurvreg"`
  and `"orsf"` learners, and refitted the model with none of the tuned
  parameters, so `compare(tune = TRUE)` silently compared default fits.
  Every built-in tuner now returns its scored grid and the refit uses
  the best configuration. A configuration that fails now gives a warning
  naming the learner and its parameter values and is left out; if every
  configuration fails,
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) stops.
  Before, failed configurations were dropped without a word, or one
  failure aborted the whole search.

- [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
  gains `method = "roc"` (time-dependent ROC curve) and `method = "dca"`
  (decision curve), both at one `eval_time`.

- [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) works
  on every
  [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
  result and returns plain column names. It failed for `"survcluster"`,
  `"roc"`, `"dca"` and `"tree_surrogate"` (“cannot coerce class … to a
  data.frame”), and for `"pdp"`, `"ice"`, `"ale"` and `"calibration"` it
  pasted the parts of the result side by side, recycling the shorter
  ones, under names such as `results.time` or
  `integrated.integrated_surv`. Now `"pdp"` and `"ice"` give one row per
  feature value, observation and time (`feature`, `feature_value`, `id`,
  `time`, `surv`, `type`), `"ale"` one row per feature value and time
  (`feature`, `feature_value`, `time`, `ale`), `"calibration"` the
  calibration table with its `eval_time`, `"survcluster"` one row per
  observation and time (`id`, `cluster`, `medoid`, `time`, `surv`),
  `"roc"` and `"dca"` the curve with its `eval_time`, and
  `"tree_surrogate"` one row per time and feature (`splits`,
  `r_squared`). For `"pdp"` and `"ale"`,
  `as.data.frame(x, what = "integrated")` returns the curve integrated
  over time (`feature`, `feature_value`, `integrated`).
  [`print()`](https://rdrr.io/r/base/print.html) shows the first rows of
  this data frame instead of the raw result list.

- [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md)
  now returns `learner` and `tunable` only, and
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md)
  returns the
  [`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md)
  method ids with a one-line summary.

- [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) no
  longer accepts engine arguments through `...`; pass them in
  `spec = list(...)`. Extra arguments now give an error instead of a
  warning.

- Fewer silent failures.
  [`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)
  warns with the number of failed bootstrap replicates and the first
  error (and stops if all fail) instead of quietly building the interval
  from the rest.
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) reports
  why it could not read a formula. The survival super learner gives the
  error when it drops a base learner. The `blackboost` and `cforest`
  learners no longer hide engine warnings; `cforest` caps `mtry` at the
  number of predictors itself, as ‘party’ would.

- The `stpm2` learner is removed: `fit_stpm2()` and `predict_stpm2()`
  are gone and ‘rstpm2’ is no longer a dependency. Use `fit_rpsurv()` /
  `fit(model = "rpsurv")`, which fits the same Royston-Parmar model on
  the hazard, odds, or probit scale and adds nonlinear effects (`nle`).

- The `orsf` tuner passed its `min_events` grid column to ‘aorsf’ as
  `n_split`, the number of candidate cut points, not the minimum number
  of events. The column is now `split_min_events` and sets that minimum.

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
  `tune_survsvm()` grid is `c(0.01, 0.1, 1)`, matched to the new loss
  scale. The fitted object’s `engine` is now `"survalis"`.

- The `rsf` learner failed with “formula is incorrectly specified” for a
  [`survival::Surv()`](https://rdrr.io/pkg/survival/man/Surv.html)
  outcome, so the survival super learner dropped it. The `rsf`, `orsf`,
  `blackboost`, `bnnsurv` and `flexsurvreg` model objects also recorded
  the wrong `time` and `status` column names.

- The `bnnsurv` learner (bagged k-nearest-neighbour Kaplan-Meier) now
  has its own engine and no longer depends on ‘bnnSurvival’, which is
  unmaintained. ‘bnnSurvival’ predicted through
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
  and forking a session whose multithreaded BLAS (e.g. MKL) is already
  running could hang `predict_bnnsurv()`. The new engine runs in the
  current process and also: one-hot encodes factors (they were integer
  codes inside the Mahalanobis distance) and matches new data to the
  training levels; uses the pseudo-inverse covariance, so constant or
  collinear columns in a bootstrap sample are handled; draws feature
  subsets from all predictors (a factor counts as one); splits weight
  equally between neighbours tied at the k-th distance, so results no
  longer depend on row order; and reads predictions at arbitrary times
  off the step function instead of interpolating linearly. With numeric
  predictors and default settings it gives the same curves as
  ‘bnnSurvival’ for the same seed, about ten times faster.
  `metric = "euclidean"` (standardized columns) is new. The fitted
  object’s `engine` is now `"survalis"`.

- survalis now declares `R (>= 4.4.0)` and `fastgbm (>= 0.7.0)`. It
  claimed R 4.1, but ‘fastgbm’ up to 0.6.1 needed R 4.5.0 and ‘rms’
  (imported through ‘pec’) needs R 4.4.0, so installation failed on
  older R.

- Fitting no longer floods the console with warnings. The `glmnet`
  learner passes `cox.ties = "breslow"`, its current default, so
  ‘glmnet’ 5.0 does not warn on every fit (this gave most of the
  warnings in a tuning or comparison run). The `rpsurv` learner muffles
  the “standard errors unavailable” warning of ‘rpsurv’ 0.8.1, since
  survalis uses only the predictions. `flexsurvreg` predictions at a
  single time no longer warn about an unknown `.pred` column.

- The `survdnn` learner is removed, and with it the ‘survdnn’ and
  ‘torch’ dependencies. ‘torch’ needs a separate LibTorch download that
  often failed to install. Use the `densemlp` learner
  (`fit(..., model = "densemlp")`), which fits deep survival models with
  a native C++ backend. `fit_survdnn()`, `predict_survdnn()`, and
  `tune_survdnn()` are removed.

- survalis no longer defines a `plot.formula()` method. It replaced
  [`graphics::plot.formula()`](https://rdrr.io/r/graphics/plot.formula.html)
  for every package as soon as survalis was loaded (“Registered S3
  method overwritten”). `plot_survcurve()` is no longer exported either:
  for a Kaplan-Meier curve use `plot(survival::survfit())` or a
  dedicated package such as ‘ggsurvfit’.
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on a
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md) result
  draws the model’s predicted survival curves.

## survalis 1.6.4

Bug fixes found while benchmarking every learner on 31 survival
datasets.

- `fit_ranger()` drops predictors that are constant in the data it is
  fitted on. With `respect.unordered.factors = "order"`, ranger stopped
  with “invalid type (NULL) for variable ‘smry\$strata’” on a factor
  with a single level, which happens in tuning folds of data with many
  rare binary predictors (`metabric`).
- `fit_glmnet()` and `predict_glmnet()` accept a single predictor.
  glmnet stops with “x should be a matrix with 2 or more columns”, so
  the predictor matrix is padded with a constant zero column.
- `fit_aalen()` standardizes numeric covariates and `predict_aalen()`
  applies the same scaling.
  [`timereg::aalen()`](https://rdrr.io/pkg/timereg/man/aalen.html)
  silently returned an all-zero fit when covariates were on very
  different scales, so the predicted survival was 1 for every subject.
- `predict_flexsurvreg()`, `predict_bnnsurv()`, and `predict_cforest()`
  accept a single evaluation time. A one-time result is now a one-column
  matrix (it was a one-row matrix, or a flat table for `flexsurvreg`),
  so the single-time metrics `brier` and `ece` can be computed for these
  learners.
- `fit_survmetalearner()` can be called with `formula` and `data` only,
  as `benchmark()` and
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
  do. It then fits the base learners (`learners`, default `coxph`,
  `glmnet`, `rsf`), builds out-of-fold predictions, and learns the
  weights from them. `predict_survmetalearner()` interpolates the
  weights when the requested times differ from the fitted grid.
- `fit_xgboost()` accepts `max_depth` and `eta`. `tune_xgboost()` passed
  them, so every configuration failed silently and
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) then
  stopped with an unrelated “Unknown columns: ibs”.
- [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) stops
  with a clear message when every configuration fails during
  cross-validation.
- `tune_survsvm()` has a default `param_grid`, and `tune_survdnn()`
  accepts its own default: an argument check rejected the default
  because [`missing()`](https://rdrr.io/r/base/missing.html) is always
  true for it.
- Errors raised inside parallel workers (`ncores > 1`) in
  [`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
  [`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
  and `cv_survlearner()` are now reported.
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
  `fit_survsvm.R`, `evaluation.R` (`cv_summary()`), `survmat-helpers.R`
  (`plot_survmat()`), and `benchmark_default_survlearners.R`
  (`summarise_benchmark()`/`summarize_benchmark_results()`/
  `best_survlearner()`) are replaced with base R and
  [`basetable::aggregate()`](https://rdrr.io/pkg/basetable/man/aggregate.html)/`count()`/`tolong()`/`towide()`/`orderrows()`/
  `applyby()`. `data.table` dropped from `Imports` entirely.
  `cv_summary()` and `score_survmodel()` etc. now return plain
  `data.frame`s instead of `data.table`s (no downstream code relied on
  data.table-specific semantics).
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
  Also fixed a latent single-bracket row-filter in `tune_bnnsurv()`
  (`results[!results$failed]`) that only worked under `data.table`’s
  `DT[i]` row-selection semantics; now explicit
  (`results[!results$failed, ]`). `data.table` remains in use elsewhere
  in the package pending further migration.

## survalis 1.6.0

- Bug fix: `survmat_to_rmst()` (and everything built on it:
  `predict(type = "rmst")`,
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
  `nbsurv` package). `fit_nbsurv()`/`predict_nbsurv()` added (no tuner).
  Registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md);
  24 learners total.

## survalis 1.4.0

- New learner: `rpsurv` (Royston-Parmar flexible parametric survival,
  via the `rpsurv` package). `fit_rpsurv()`/`predict_rpsurv()` added (no
  tuner), restricted cubic splines with optional time-varying covariate
  effects (`tve`) and a choice of `scale` (hazard/odds/normal).
  Registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md);
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
- `default_times()`: a recorded, data-dependent evaluation grid used
  whenever `times` is omitted, so runs stay reproducible.
- Resampling constructors
  [`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md)
  (with `repeats`),
  [`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
  [`group_cv()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
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
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md);
  22 learners total.

### Deprecations

- The granular functions (`fit_*()`, `predict_*()`, `tune_*()`,
  `cv_survlearner()`, `score_survmodel()`, `benchmark*()`,
  `summarise_benchmark()`, `estimate_rmst()`, …) remain exported and
  fully supported. They are now marked internal in the manual and emit a
  one-time per-session note. See `?\`survalis-deprecated\`\`. They are
  scheduled to become unexported in a future release.

## survalis 1.2.0

- Added `compute_survcluster()` and `plot_survcluster()`: unsupervised
  clustering of predicted individualized survival curves via
  [`unsurv::unsurv()`](https://rdrr.io/pkg/unsurv/man/unsurv.html) (PAM
  on a weighted curve representation), registered in
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md).
  Added `unsurv` to `Imports` and its paper’s citation (El Badisy 2026,
  `<doi:10.1093/bioadv/vbag218>`) to `Description`.

## survalis 1.1.0

- Added `screen_fdr()`: univariate-Cox variable screening with
  Benjamini-Hochberg FDR control, with a `minscreen` backfill guarantee.
- Added `estimate_rmst()` and `plot_rmst()`: standardized RMST
  estimation via g-computation, with percentile bootstrap inference.
  Dispatch between a marginal estimate and a two-arm causal-style
  contrast is by the `trt_col` argument (`NULL` = marginal), matching
  `plot_survmat()`’s existing `group` dispatch convention rather than
  separate function names.
- Added `compute_shap_matrix()` and `plot_shap_beeswarm()`: a beeswarm
  view of SHAP contributions time-integrated over the full prediction
  horizon, built on `compute_shap()`’s existing per-time Shapley
  decomposition. `compute_shap()` gained an explicit `features` argument
  (defaulting to
  `setdiff(colnames(newdata), c(model$time, model$status))`), replacing
  a hardcoded `c("time", "status", "event")` exclusion, for consistency
  with `compute_interactions()`’s existing `features` argument.
- Added `roc_survmat()` (full time-dependent ROC curve, IPCW-weighted
  like `auc_survmat()`) and `dca_survmat()` (decision curve analysis /
  net benefit for right-censored outcomes, using the same `Ĝ(t)` IPCW
  estimator as `brier()` and `auc_survmat()`), each with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)/`plot_roc()`/`plot_dca()`
  methods.

## survalis 1.0.0

CRAN release: 2026-08-22

- Added `CoxBoost` as a new survival learner: `fit_coxboost()`,
  `predict_coxboost()`, and `tune_coxboost()`, wrapping the `CoxBoost`
  package’s likelihood-based componentwise boosting for the Cox partial
  likelihood. Registered in
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md)
  and usable directly through `benchmark()`.
- Wrapped `tune_glmnet()`’s example in `\donttest{}` (per CRAN
  win-builder pretest feedback: it exceeded the 5s/10s example runtime
  limit on Debian/Windows).

## survalis 0.9.0

- Added a `title` argument to every plotting function (`plot_ale()`,
  `plot_pdp()`, `plot_shap()`, `plot_interactions()`,
  `plot_calibration()`, `plot_counterfactual()`, `plot_surrogate()`,
  `plot_tree_surrogate()`, `plot_varimp()`, `plot_survmat()`,
  `plot_survmetalearner_weights()`, `plot_benchmark()`, `cv_plot()`,
  `plot_survcurve()`). If omitted, the previous automatically generated
  title is used (fully backward compatible); pass `title = NULL` to omit
  the title entirely, for journals that require caption-only figures.

## survalis 0.8.9

- Added a dashed zero-reference line to every plot that can show both
  negative and positive values (previously only `plot_shap()`’s bar
  chart had one): `plot_ale()` (both `per_time` and `integrated`),
  `plot_shap()`’s time-dependent curve, `plot_counterfactual()`,
  `plot_surrogate()`, and `plot_varimp()`.

## survalis 0.8.8

- `plot_varimp()` now draws a boxplot of the full per-repetition
  permutation distribution for each feature instead of a single point
  estimate with no uncertainty shown. `compute_varimp()` retains the
  per-repetition values (previously discarded after computing summary
  stats) as a `"raw_scores"` attribute on its return value, used by
  `plot_varimp()`; falls back to the previous point-plot behavior for
  objects without that attribute (e.g. hand-built summary tables).

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
  ggplot2’s default hue palette everywhere except `plot_survcurve()`;
  they now all use `scale_color_survalis()`/`scale_fill_survalis()`.
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

- `summarise_benchmark()` now rounds `mean`/`sd`/`se`/`lower`/`upper` to
  3 decimals by default (`digits = 3` argument), matching
  `cv_summary()`’s convention; it previously returned full-precision
  values (e.g. `sd = 0.014142136`).

## survalis 0.8.3

- Normalized `base_size` to the shared `theme_survalis()` default (13)
  everywhere: `compute_shap.R`, `fit_survmetalearner.R`’s weights plot,
  and `plot_benchmark()` previously carried over a hardcoded
  `base_size = 14` from before the theme retrofit (PR
  [\#22](https://github.com/ielbadisy/survalis/issues/22)).

## survalis 0.8.2

- Migrated `compute_shap()`, `compute_shap_mean()`,
  `compute_calibration()`, and `plot_survmetalearner_weights()` from
  dplyr/tidyr to data.table. This was the remaining unqualified
  (non-`dplyr::`-prefixed) dplyr/tidyr usage missed by the original
  file-by-file migration inventory.
- `dplyr`, `tidyr`, `purrr`, and `tibble` are now fully removed from
  `Imports`: the package no longer depends on any tidyverse package for
  its data-manipulation code; `data.table` is the sole engine
  throughout.

## survalis 0.8.1

- Migrated the nested/tuned-CV half of
  `benchmark_default_survlearners.R`
  (`benchmark_default_survlearners()`, `.nested_surv_*` helpers,
  `summarise_benchmark()`, `summarize_benchmark_results()`,
  `best_survlearner()`) from dplyr/tidyr/tibble to data.table. This
  completes the file-by-file dplyr/tidyverse -\> data.table migration
  across the package. Reuses the shared `.score_metrics()` helper in
  `.nested_surv_score_predictions()` instead of a fifth duplicated copy
  of the per-metric scoring
  [`switch()`](https://rdrr.io/r/base/switch.html) block.

## survalis 0.8.0

- Migrated `cv_survmetalearner()` from dplyr/tidyr/purrr/tibble to
  data.table, reusing the shared `.score_metrics()` helper instead of a
  third duplicated copy of the per-metric scoring
  [`switch()`](https://rdrr.io/r/base/switch.html) block.

## survalis 0.7.27

- Migrated
  [`list_survlearners()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  [`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md),
  and
  [`list_metrics()`](https://elbadisyimad.com/survalis/reference/survalis-deprecated.md)
  from
  [`tibble::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
  to
  [`data.table::data.table()`](https://rdrr.io/pkg/data.table/man/data.table.html)
  (cosmetic descriptive tables, no behavior change).

## survalis 0.7.26

- Migrated `plot_survmat()` from dplyr/tidyr to data.table
  ([`data.table::melt()`](https://rdrr.io/pkg/data.table/man/melt.data.table.html)
  for the wide-to-long reshape, `[.data.table` with `by =` for the
  group/time summaries).

## survalis 0.7.25

- Migrated `compute_varimp()` from dplyr/purrr/tibble to data.table.

## survalis 0.7.24

- Migrated `tune_survdnn()` from dplyr/tidyr/purrr/tibble to data.table,
  including the list-valued `hidden` architecture parameter (grid
  expansion built manually over unique levels per parameter, since
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html)
  doesn’t support list-valued columns).

## survalis 0.7.23

- Migrated `tune_survsvm()` from dplyr/tidyr/purrr/tibble to data.table
  (grid expansion via
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html)); the
  CRAN-motivated refit-fallback behavior (PR
  [\#14](https://github.com/ielbadisy/survalis/issues/14)) is unchanged
  and still verified passing.

## survalis 0.7.22

- Migrated `tune_glmnet()` from dplyr/tidyr/purrr/tibble to data.table,
  including its
  [`tidyr::crossing()`](https://tidyr.tidyverse.org/reference/expand.html)
  grid-expansion step (now
  [`data.table::CJ()`](https://rdrr.io/pkg/data.table/man/J.html), a
  sorted-unique cross join with matching semantics).

## survalis 0.7.21

- Migrated `tune_selectcox()` from dplyr/tidyr/purrr/tibble to
  data.table.

## survalis 0.7.20

- Migrated `tune_flexsurvreg()` from dplyr/tidyr/purrr/tibble to
  data.table. Added `.map_rbind_dt()` in `R/dt-utils.R`, a data.table
  replacement for
  [`purrr::map_dfr()`](https://purrr.tidyverse.org/reference/map_dfr.html)
  over a plain vector (as opposed to `.pmap_rbind_dt()` for a parameter
  grid).

## survalis 0.7.19

- Migrated `tune_blackboost()` from dplyr/tidyr/purrr/tibble to
  data.table.

## survalis 0.7.18

- Migrated `tune_bart()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.17

- Migrated `tune_rsf()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.16

- Migrated `tune_orsf()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.15

- Migrated `tune_cforest()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.14

- Migrated `tune_bnnsurv()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.13

- Migrated `tune_xgboost()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.12

- Migrated `tune_ranger()` from dplyr/tidyr/purrr/tibble to data.table.

## survalis 0.7.11

- Fixed `benchmark_tuned_survlearners()` silently dropping learners
  whose `tune_*()` has been migrated to data.table:
  `tuning_results[1, cols, drop = FALSE]` relies on data.frame `[`
  semantics, but data.table’s `[` does not select columns when `cols` is
  a variable (only a literal character vector triggers that); it was
  returning the raw column-name vector instead of a one-row parameter
  table, causing every fold to error and the learner to be dropped with
  a warning. Added `.select_cols()` in `R/dt-utils.R` that selects
  columns correctly for either a data.frame/tibble or a data.table, so
  the nested-tuning code stays agnostic to migration state.

## survalis 0.7.10

- Retrofitted all existing plotting functions to use `theme_survalis()`
  instead of ad-hoc `theme_minimal()` calls with inconsistent
  `base_size` values, so package figures share one consistent visual
  style: `plot_ale()`, survmat plotting helpers,
  `plot_surrogate()`/`plot_shap()`, `cv_plot()`, `plot_pdp()`,
  `plot_interactions()`, `plot_counterfactual()`,
  `plot_tree_surrogate()`, `plot_benchmark()`, `plot_varimp()`,
  `fit_survmetalearner()`’s plotting path. Hardcoded per-plot colors are
  unchanged in this pass (follow-up).

## survalis 0.7.9

- Added `plot_survcurve()`, a survminer-inspired Kaplan-Meier curve with
  confidence ribbon, log-rank p-value annotation, and an aligned
  number-at-risk table, built natively on `theme_survalis()` /
  `scale_color_survalis()` (no dependency on the survminer package). The
  risk table uses `patchwork` (added to Suggests).

## survalis 0.7.8

- Added `theme_survalis()`, a shared ggplot2 theme, and
  `scale_color_survalis()`/`scale_fill_survalis()`, colorblind-friendly
  discrete scales (ColorBrewer “Dark2”-based, recycled/interpolated
  beyond 8 levels), for a consistent visual style across package
  figures. Not yet applied to existing plotting functions (follow-up).

## survalis 0.7.7

- Migrated `tune_rpart()` from dplyr/tidyr/purrr/tibble to data.table.
  Added internal data.table helpers (`R/dt-utils.R`) shared by the
  ongoing file-by-file migration of the remaining `tune_*()`/`fit_*()`
  grid-search code (tracked in the project TODO).

## survalis 0.7.6

- Added `benchmark()`, the single entry point for comparing survival
  learners: `tune = FALSE` (default) dispatches to
  `benchmark_default_survlearners()`, `tune = TRUE` dispatches to
  `benchmark_tuned_survlearners()` (nested CV with per-learner tuning).
  `benchmark_default_survlearners()` and
  `benchmark_tuned_survlearners()` remain available directly.

## survalis 0.7.5

- Migrated the core CV/metrics engine (`cv_survlearner()`,
  `cv_summary()`, `score_survmodel()`) from dplyr/tidyr/tibble/purrr to
  data.table; `data.table` is now the main data-manipulation engine
  going forward, replacing the dplyr/tidyverse approach used previously.
  This is a first step; the remaining tuning/fitting code
  (`tune_*()`/`fit_*()`) still uses dplyr/purrr and will be migrated
  incrementally.
- `cv_summary()` and `score_survmodel()`/`cv_survlearner()` metric
  values now round to 3 decimals by default (`digits = 3` argument on
  `cv_summary()`).
- These functions now return `data.table` objects instead of tibbles
  ([`is.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) still
  holds; code relying on `tibble`-specific behavior should call
  [`tibble::as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
  explicitly).

## survalis 0.7.4

- Added `timeroc_survmat()`, a vectorized cumulative/dynamic
  time-dependent AUC curve over a vector of evaluation times (Uno et
  al. 2007 / Heagerty and Zheng 2005 estimator), matching
  `timeROC::timeROC(weighting = "marginal")` to ~1e-3.
- Fixed `auc_survmat()` to define cases as events strictly before
  `t_star` (`time < t_star`) rather than `time <= t_star`, matching the
  canonical Uno/timeROC definition; this corrects the `"auc"` metric
  used throughout `score_survmodel()`, `benchmark_*()`, and
  `fit_survmetalearner()`.
- Added `timeROC` to Suggests for numerical validation in tests.

## survalis 0.7.3

- Fixed `plot_pdp()` per-time facets using free y-axis scales, which
  made survival probability panels visually incomparable across facets;
  now uses a fixed `[0, 1]` scale via `coord_cartesian()`.

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
- Added `plot_counterfactual()` for visualizing counterfactual
  recommendations.
- Added `plot_survmat()` for plotting predicted survival curves,
  including grouped summaries.
- Hardened documentation and examples for CRAN checks.
- Updated dependency declarations so exported learners and CV/tuning
  helpers work after installation.
- Added clearer runtime checks for `survdnn` when the LibTorch runtime
  is missing.
