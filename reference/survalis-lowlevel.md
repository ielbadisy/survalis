# Low-level functions in survalis

survalis 1.3.0 introduced the seven-verb interface –
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md),
[`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md)
– plus the resampling constructors
([`cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`holdout()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`group_cv()`](https://elbadisyimad.com/survalis/reference/resampling.md),
[`bootstrap()`](https://elbadisyimad.com/survalis/reference/resampling.md)).

## Details

The granular functions listed below are the machinery those verbs are
built on. They stay exported and supported for direct use, for example
in simulation code that calls one learner many times.

## Mapping

|  |  |
|----|----|
| **Low-level** | **Verb** |
| `fit_<learner>()` | `fit(formula, data, model = "<learner>")` |
| `predict_<learner>()` | `predict(fit, newdata, times, type)` |
| `tune_<learner>()` | `tune(formula, data, model = "<learner>", grid)` |
| [`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md), [`cv_summary()`](https://elbadisyimad.com/survalis/reference/cv_summary.md), [`cv_plot()`](https://elbadisyimad.com/survalis/reference/cv_plot.md) | `evaluate(...)` |
| [`score_survmodel()`](https://elbadisyimad.com/survalis/reference/score_survmodel.md) | `evaluate(fit, resampling = ...)` |
| [`benchmark()`](https://elbadisyimad.com/survalis/reference/benchmark.md), [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md), [`benchmark_tuned_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_tuned_survlearners.md) | `compare(...)` |
| [`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md), [`summarize_benchmark_results()`](https://elbadisyimad.com/survalis/reference/summarize_benchmark_results.md), [`best_survlearner()`](https://elbadisyimad.com/survalis/reference/best_survlearner.md) | [`summary()`](https://rdrr.io/r/base/summary.html) / [`print()`](https://rdrr.io/r/base/print.html) on a `survalis_compare` |
| `compute_*()` / `plot_*()` | `interpret(fit, method = ...)` then [`plot()`](https://rdrr.io/r/graphics/plot.default.html) |
| [`estimate_rmst()`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md) | `estimate(..., estimand = "RMST_diff")` |
