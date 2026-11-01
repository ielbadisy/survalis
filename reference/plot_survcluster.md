# Plot Clustered Survival Curves

Plots the individualized survival curves from
[`compute_survcluster()`](https://elbadisyimad.com/survalis/reference/compute_survcluster.md),
colored by cluster membership, with medoid curves overlaid.

## Usage

``` r
plot_survcluster(survcluster_output, title)
```

## Arguments

- survcluster_output:

  Output list returned by
  [`compute_survcluster()`](https://elbadisyimad.com/survalis/reference/compute_survcluster.md).

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A `ggplot2` object showing per-subject survival curves colored by
cluster, with cluster medoid curves overlaid.

## See also

[`compute_survcluster()`](https://elbadisyimad.com/survalis/reference/compute_survcluster.md)

## Examples

``` r
if (requireNamespace("unsurv", quietly = TRUE) &&
    requireNamespace("cluster", quietly = TRUE)) {
  mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
  time_grid <- seq(10, 300, length.out = 20)
  cl <- compute_survcluster(mod, newdata = veteran, times = time_grid, K = 2, seed = 1)
  plot_survcluster(cl)
}

```
