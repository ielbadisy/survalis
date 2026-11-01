# Plot Benchmark Distributions Across Learners

Produces box‑and‑jitter plots of CV metric values per learner, faceted
by metric for quick visual comparison.

## Usage

``` r
plot_benchmark(benchmark_results, title)
```

## Arguments

- benchmark_results:

  A data frame from
  [`benchmark_default_survlearners()`](https://elbadisyimad.com/survalis/reference/benchmark_default_survlearners.md)
  with columns `learner`, `metric`, and `value`.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object.

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
plot_benchmark(res)

```
