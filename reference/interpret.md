# Interpret a fitted survival learner

A single entry point to the package's model-agnostic explanation and
diagnostic methods.
[`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md)
lists the `method` ids with a one-line description of each.

## Usage

``` r
interpret(fit, method, newdata = NULL, data = NULL, times = NULL, ...)
```

## Arguments

- fit:

  A `survalis_fit` (from
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md)), a
  `survalis_tune`, or a legacy `mlsurv_model`.

- method:

  One of `"pdp"`, `"ice"`, `"ale"`, `"shap"`, `"shap_matrix"`,
  `"varimp"` / `"permute"`, `"interaction"`, `"counterfactual"`,
  `"surrogate"`, `"tree_surrogate"`, `"calibration"`, `"roc"`, `"dca"`,
  `"survcluster"`.

- newdata:

  Instances to explain (SHAP, ALE, counterfactual, surrogate,
  survcluster) or to assess (`roc`, `dca`; these need the outcome
  columns). Defaults to `data`, so `roc` and `dca` are apparent unless
  you pass held-out data.

- data:

  Background data for the method. Defaults to the training data.

- times:

  Evaluation grid; `NULL` picks one from the data quantiles. Methods
  that need a single time (`calibration`, `roc`, `dca`) take `eval_time`
  in `...` or the sole `times` value.

- ...:

  Method-specific arguments, e.g. `feature` (pdp / ice / ale),
  `target_time` (counterfactual, surrogate), `eval_time` (calibration,
  roc, dca), `thresholds` (dca), `K` (survcluster). For `"shap"`, the
  aggregation over time (`"meanabs"` or `"integral"`) is passed as
  `shap_method`.

## Value

A `survalis_interpret` object holding the method's `result`, with
`method`, `times`, and a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the result frame. For `"survcluster"` it is in long format, one row per
observation and time: `id` (row of `newdata`), `cluster`, `medoid`
(`TRUE` for the cluster's medoid curve), `time` and predicted survival
`surv`. For `"roc"` and `"dca"` it is the curve with an `eval_time`
column (and `auc` for `"roc"`); for `"tree_surrogate"` one row per time
and feature with the number of `splits` on it and the tree's
`r_squared`. For `"pdp"` / `"ice"` it is one row per curve point
(`feature`, `feature_value`, `id` (`NA` on PDP rows), `time`, `surv`,
`type`), for `"ale"` one row per feature value and time (`feature`,
`feature_value`, `time`, `ale`), and for `"calibration"` the binned
table with its `eval_time`. For `"pdp"` and `"ale"`,
`as.data.frame(x, what = "integrated")` returns the curve integrated
over time instead (`feature`, `feature_value`, `integrated`).

## See also

For scalar summaries use
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md):
`metrics = "ece"` summarises calibration and `metrics = "auc"`
summarises the ROC curve.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
         spec = list(num.trees = 100))
interpret(m, "pdp", feature = "karno", times = c(90, 180))
#> <survalis_interpret> method = pdp  |  model = ranger
#>   feature: karno 
#>   grid: 2 times in [90, 180]
#>   $result:
#>   feature feature_value id time      surv type
#> 1   karno      10.00000  1   90 0.1884062  ice
#> 2   karno      14.68421  1   90 0.1884062  ice
#> 3   karno      19.36842  1   90 0.1875603  ice
#> 4   karno      24.05263  1   90 0.1875603  ice
#> 5   karno      28.73684  1   90 0.2404616  ice
interpret(m, "varimp", times = c(90, 180))
#> <survalis_interpret> method = varimp  |  model = ranger
#>   grid: 2 times in [90, 180]
#>   $result:
#>    feature importance importance_05 importance_95 scaled_importance
#> 1    karno     0.1408       0.11545       0.17390         100.00000
#> 2 celltype     0.0672       0.05490       0.08305          47.72727
#> 3      age     0.0439       0.03280       0.05395          31.17898
#> 4      trt     0.0000       0.00000       0.00000           0.00000
#> 5 diagtime     0.0000       0.00000       0.00000           0.00000
plot(interpret(m, "roc", eval_time = 180))
```
