# Select the Best Survival Learner by a Given Metric

Extracts the top‑performing learner(s) under a chosen metric from
benchmark results, using the average value across folds.

## Usage

``` r
best_survlearner(benchmark_results, metric, maximize = NULL)
```

## Arguments

- benchmark_results:

  A data frame with columns `learner`, `metric`, and `value`.

- metric:

  Character name of the metric to optimize (e.g., `"cindex"`, `"ibs"`).
  Must exist in `benchmark_results$metric`.

- maximize:

  Logical; whether to maximize the metric. If `NULL` (default), the
  function maximizes for concordance‑like metrics (`"cindex"`) and
  minimizes for error‑like metrics (`"ibs"`, `"brier"`, `"iae"`,
  `"ise"`, `"ece"`).

## Value

A tibble with columns `learner`, `metric`, and the selected average
`value` for the best learner(s). Ties are returned as multiple rows.

## See also

[`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md),
[`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md)

## Examples

``` r
res <- data.frame(
  learner = c("coxph", "coxph", "rpart", "rpart"),
  metric = c("cindex", "ibs", "cindex", "ibs"),
  value = c(0.64, 0.19, 0.60, 0.23)
)
best_survlearner(res, metric = "cindex")
#> # basetable: 1 x 3
#>   learner metric value
#> 1   coxph cindex  0.64
best_survlearner(res, metric = "ibs")
#> # basetable: 1 x 3
#>   learner metric value
#> 1   coxph    ibs  0.19
```
