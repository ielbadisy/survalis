# Summarise Benchmark Results (Mean SD with Wald CI)

Aggregates cross‑validated benchmark results by learner and metric,
reporting mean, standard deviation, standard error, and an approximate
95% Wald confidence interval.

## Usage

``` r
summarise_benchmark(benchmark_results, digits = 3)
```

## Arguments

- benchmark_results:

  A data frame produced by
  [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md),
  containing at least `learner`, `metric`, and `value`.

- digits:

  Integer number of decimal places for `mean`, `sd`, `se`, `lower`,
  `upper` (default `3`).

## Value

A data.frame with columns `learner`, `metric`, `mean`, `sd`, `n`, `se`,
`lower`, `upper`.

## See also

[`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md),
[`plot_benchmark()`](https://elbadisyimad.com/survalis/reference/plot_benchmark.md),
[`summarize_benchmark_results()`](https://elbadisyimad.com/survalis/reference/summarize_benchmark_results.md)

## Examples

``` r
res <- data.frame(
  learner = c("coxph", "coxph", "rpart", "rpart"),
  metric = c("cindex", "ibs", "cindex", "ibs"),
  value = c(0.64, 0.19, 0.60, 0.23)
)
summarise_benchmark(res)
#>   learner metric mean sd n se lower upper
#> 1   coxph cindex 0.64 NA 1 NA    NA    NA
#> 2   coxph    ibs 0.19 NA 1 NA    NA    NA
#> 3   rpart cindex 0.60 NA 1 NA    NA    NA
#> 4   rpart    ibs 0.23 NA 1 NA    NA    NA
```
