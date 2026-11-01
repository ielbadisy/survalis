# Plot ALE Curves for Survival Models

Visualizes ALE results produced by
[`compute_ale()`](https://elbadisyimad.com/survalis/reference/compute_ale.md)
either as per-time curves (one curve per evaluation time) or as an
integrated curve averaged across times.

## Usage

``` r
plot_ale(
  ale_result,
  feature,
  which = c("per_time", "integrated"),
  smooth = FALSE,
  title
)
```

## Arguments

- ale_result:

  A list returned by
  [`compute_ale()`](https://elbadisyimad.com/survalis/reference/compute_ale.md).

- feature:

  Character name of the feature (for axis labeling only).

- which:

  Either `"per_time"` to plot one ALE curve per time, or `"integrated"`
  to plot the time-averaged ALE (requires multiple `times`).

- smooth:

  Logical; if `TRUE`, overlays a LOESS smooth on the plotted curve(s).

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A `ggplot2` object.

## Details

Per-time plots show how the feature's local effect varies across
different evaluation times. The integrated plot summarizes the average
effect over the supplied time grid (simple mean across times of the
centered ALE values).

## See also

[`compute_ale()`](https://elbadisyimad.com/survalis/reference/compute_ale.md),
[`compute_pdp()`](https://elbadisyimad.com/survalis/reference/compute_pdp.md),
[`plot_pdp()`](https://elbadisyimad.com/survalis/reference/plot_pdp.md)

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
ale_res <- compute_ale(
  model = mod,
  newdata = veteran,
  feature = "karno",
  times = c(80, 160),
  grid.size = 8
)
plot_ale(ale_res, feature = "karno", which = "per_time")

plot_ale(ale_res, feature = "karno", which = "integrated", smooth = TRUE)
#> `geom_smooth()` using formula = 'y ~ x'

```
