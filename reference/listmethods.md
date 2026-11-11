# List the interpretability methods

The ids accepted by the `method` argument of
[`interpret()`](https://elbadisyimad.com/survalis/reference/interpret.md).

## Usage

``` r
listmethods()
```

## Value

A data frame with columns `method` and `summary`.

## Examples

``` r
listmethods()
#>            method
#> 1             pdp
#> 2             ice
#> 3             ale
#> 4            shap
#> 5     shap_matrix
#> 6          varimp
#> 7         permute
#> 8     interaction
#> 9  counterfactual
#> 10      surrogate
#> 11 tree_surrogate
#> 12    calibration
#> 13            roc
#> 14            dca
#> 15    survcluster
#>                                                                    summary
#> 1      Partial dependence of predicted survival on one feature, over time.
#> 2               Individual conditional expectation curves for one feature.
#> 3                               Accumulated local effects for one feature.
#> 4                   Shapley values for one instance, aggregated over time.
#> 5                       Shapley values for many instances (beeswarm plot).
#> 6                   Permutation variable importance on the survival scale.
#> 7                                                         Alias of varimp.
#> 8                   Friedman H-statistic interaction strength per feature.
#> 9  Smallest feature changes that move predicted survival at a target time.
#> 10               Local linear surrogate for one instance at a target time.
#> 11                      Global decision-tree surrogate of the predictions.
#> 12                                          Calibration curve at one time.
#> 13                                   Time-dependent ROC curve at one time.
#> 14          Decision curve (net benefit over risk thresholds) at one time.
#> 15            Clusters of subjects with similar predicted survival curves.
```
