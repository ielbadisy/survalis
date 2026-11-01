# List the evaluation metrics

The ids accepted by the `metrics` argument of
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md)
and
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md)
and the `metric` argument of
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md).

## Usage

``` r
list_metrics()
```

## Value

A data.frame with one row per metric.

## Details

Columns:

- `metric`: short metric name used throughout the package.

- `direction`: whether higher or lower values are better.

- `summary`: brief description.

- `range`: typical value range.

## Examples

``` r
list_metrics()
#>   metric direction
#> 1 cindex  maximize
#> 2    auc  maximize
#> 3  brier  minimize
#> 4    ibs  minimize
#> 5    iae  minimize
#> 6    ise  minimize
#> 7    ece  minimize
#>                                                                    summary
#> 1                Harrell-style concordance index for survival predictions.
#> 2     Cumulative/dynamic time-dependent AUC at a selected evaluation time.
#> 3 Brier Score at specified evaluation time(s) (IPCW-weighted when needed).
#> 4     Integrated Brier Score over an evaluation time grid (IPCW-weighted).
#> 5                Integrated absolute error against the Kaplan-Meier curve.
#> 6                 Integrated squared error against the Kaplan-Meier curve.
#> 7                  Expected calibration error at a single evaluation time.
#>                        range
#> 1  [0, 1] (higher is better)
#> 2  [0, 1] (higher is better)
#> 3   [0, 1] (lower is better)
#> 4   [0, 1] (lower is better)
#> 5 [0, Inf) (lower is better)
#> 6 [0, Inf) (lower is better)
#> 7   [0, 1] (lower is better)
```
