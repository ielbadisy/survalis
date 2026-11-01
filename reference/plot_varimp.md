# Plot Permutation Variable Importance

Plots the distribution of per-repetition permutation importance values
as a boxplot per feature, using either the scaled importance (default)
or the raw importance column. Falls back to a single point per feature
(with no distribution shown) for `varimp_df` objects that lack the
`"raw_scores"` attribute (e.g., hand-built summary tables).

## Usage

``` r
plot_varimp(varimp_df, use_scaled = TRUE, title)
```

## Arguments

- varimp_df:

  A data frame as returned by
  [`compute_varimp`](https://elbadisyimad.com/survalis/reference/compute_varimp.md).

- use_scaled:

  Logical; if `TRUE` (default), plot `scaled_importance` (percent). If
  unavailable, falls back to raw `importance` with a warning.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object.

## Examples

``` r
mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
imp <- compute_varimp(
  model = mod,
  times = 80,
  metric = "brier",
  n_repetitions = 3,
  seed = 1,
  subset = 40
)
plot_varimp(imp, use_scaled = TRUE)

plot_varimp(imp, use_scaled = FALSE)
```
