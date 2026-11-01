# Interpret a fitted survival learner

A single entry point to the package's model-agnostic explanation
methods. The available `method` ids are derived from
[`list_interpretability_methods()`](https://elbadisyimad.com/survalis/reference/list_interpretability_methods.md)
(plus `"ice"`, a variant of `"pdp"`, and `"permute"`, an alias of
`"varimp"`).

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
  `"surrogate"`, `"tree_surrogate"`, `"calibration"`, `"survcluster"`.

- newdata:

  Instances to explain (SHAP, ALE, counterfactual, surrogate,
  survcluster). Defaults to `data`.

- data:

  Background data for the method. Defaults to the training data.

- times:

  Evaluation grid; `NULL` uses
  [`default_times()`](https://elbadisyimad.com/survalis/reference/default_times.md).
  Methods that need a single time (`calibration`) take `eval_time` in
  `...` or the sole `times` value.

- ...:

  Method-specific arguments, e.g. `feature` (pdp / ice / ale),
  `target_time` (counterfactual), `eval_time` (calibration), `K`
  (survcluster). For `"shap"`,
  [`compute_shap()`](https://elbadisyimad.com/survalis/reference/compute_shap.md)'s
  own aggregation `method` argument (`"meanabs"` / `"integral"`) is
  passed as `shap_method`, since `method` here already names the
  `interpret()` dispatch argument.

## Value

A `survalis_interpret` object wrapping the underlying `compute_*`
result, with `method`, `times`, and a
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method that
dispatches to the paired `plot_*` helper.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the raw result frame.

## See also

For the scalar calibration summary use `evaluate(metrics = "ece")` /
[`ece_survmat()`](https://elbadisyimad.com/survalis/reference/ece_survmat.md);
`method = "calibration"` here is the calibration *curve*.

## Examples

``` r
m <- fit(Surv(time, status) ~ age + karno + celltype, veteran, "ranger",
         spec = list(num.trees = 100))
interpret(m, "pdp", feature = "karno", times = c(90, 180))
#> <survalis_interpret> method = pdp  |  model = ranger
#>   feature: karno 
#>   grid: 2 times in [90, 180]
#>   $result:
#>       results.karno results..id results.time results.surv_prob results.type
#> ice.1      10.00000           1           90         0.2069254          ice
#> ice.2      14.68421           1           90         0.2069254          ice
#> ice.3      19.36842           1           90         0.2062368          ice
#> ice.4      24.05263           1           90         0.2062368          ice
#> ice.5      28.73684           1           90         0.2554830          ice
#>       pdp_integrated.karno pdp_integrated.integrated_surv
#> ice.1             10.00000                       16.97547
#> ice.2             14.68421                       16.97547
#> ice.3             19.36842                       16.77698
#> ice.4             24.05263                       16.77698
#> ice.5             28.73684                       17.34241
interpret(m, "varimp", times = c(90, 180))
#> <survalis_interpret> method = varimp  |  model = ranger
#>   grid: 2 times in [90, 180]
#>   $result:
#>    feature importance importance_05 importance_95 scaled_importance
#> 1    karno     0.1226       0.10825       0.14000         100.00000
#> 2 celltype     0.0630       0.04650       0.07310          51.38662
#> 3      age     0.0451       0.03405       0.05465          36.78630
#> 4      trt     0.0000       0.00000       0.00000           0.00000
#> 5 diagtime     0.0000       0.00000       0.00000           0.00000
```
