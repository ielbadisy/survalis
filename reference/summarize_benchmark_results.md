# Compact Table of Mean SD by Learner and Metric

Creates a wide table summarizing each learner's performance as `mean sd`
per metric, suitable for reporting.

## Usage

``` r
summarize_benchmark_results(results, digits = 3)
```

## Arguments

- results:

  A data frame with columns `learner`, `metric`, and `value` (e.g.,
  output of
  [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md)).

- digits:

  Integer number of decimal places in the formatted summary (default
  `3`).

## Value

A wide tibble with one row per learner and one column per metric,
containing formatted strings `"mean sd"`.

## See also

[`summarise_benchmark()`](https://elbadisyimad.com/survalis/reference/summarise_benchmark.md),
[`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md)

## Examples

``` r
res <- data.frame(
  learner = c("coxph", "coxph", "rpart", "rpart"),
  metric = c("cindex", "ibs", "cindex", "ibs"),
  value = c(0.64, 0.19, 0.60, 0.23)
)
summarize_benchmark_results(res, digits = 2)
#> # basetable: 2 x 3
#>   learner   cindex      ibs
#> 1   coxph 0.64  NA 0.19  NA
#> 2   rpart 0.60  NA 0.23  NA
```
