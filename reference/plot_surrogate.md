# Plot Local Surrogate Explanation

Visualizes the signed local effects from
[`compute_surrogate`](https://elbadisyimad.com/survalis/reference/compute_surrogate.md)
as a horizontal bar chart, optionally limiting to the top `n`
contributors.

## Usage

``` r
plot_surrogate(surrogate_df, top_n = NULL, title)
```

## Arguments

- surrogate_df:

  A data frame returned by
  [`compute_surrogate()`](https://elbadisyimad.com/survalis/reference/compute_surrogate.md).

- top_n:

  Optional integer; if provided, display only the top `top_n` features
  by absolute effect.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object showing feature contributions \\\beta_j \cdot x_j\\
(positive/negative) at the target time.

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
local_expl <- compute_surrogate(
  model = mod,
  newdata = veteran[2, , drop = FALSE],
  baseline_data = veteran,
  times = c(80, 160),
  target_time = 80,
  k = 3
)
plot_surrogate(local_expl, top_n = 3)
```
