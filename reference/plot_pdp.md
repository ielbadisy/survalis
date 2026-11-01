# Plot PDP/ICE Curves for Survival Models

Plots partial dependence (PDP) and/or individual conditional expectation
(ICE) results returned by
[`compute_pdp`](https://elbadisyimad.com/survalis/reference/compute_pdp.md)
either per evaluation time or as an integrated PDP over time.

## Usage

``` r
plot_pdp(
  pdp_ice_output,
  feature,
  method = "pdp+ice",
  ids = NULL,
  which = c("per_time", "integrated"),
  alpha_ice = 0.2,
  smooth = FALSE,
  title
)
```

## Arguments

- pdp_ice_output:

  The list returned by
  [`compute_pdp()`](https://elbadisyimad.com/survalis/reference/compute_pdp.md).

- feature:

  Character scalar; the same feature analyzed in
  [`compute_pdp()`](https://elbadisyimad.com/survalis/reference/compute_pdp.md).

- method:

  One of `"pdp"`, `"ice"`, or `"pdp+ice"` indicating which curves to
  draw in per-time plots (ignored for integrated plots).

- ids:

  Optional vector of row ids (`.id`) to subset ICE curves for clarity in
  per-time plots.

- which:

  One of `"per_time"` (facet by time) or `"integrated"` (plot the
  time-integrated PDP if available).

- alpha_ice:

  Alpha transparency for ICE lines/boxes in per-time plots (default
  0.2).

- smooth:

  Logical; if `TRUE` and the feature is numeric, apply a smooth curve
  (`loess`) for integrated PDP; otherwise draw a line.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely (e.g., for journals requiring
  caption-only figures).

## Value

A ggplot2 object.

## Details

**Per-time:** For numeric features, draws ICE lines and PDP overlays per
time. For categorical features, shows ICE as boxplots per level and PDP
as point summaries. **Integrated:** Plots the PDP integrated across time
(if provided by
[`compute_pdp()`](https://elbadisyimad.com/survalis/reference/compute_pdp.md));
numeric features can be smoothed with `smooth=TRUE`.

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
pdp_age <- compute_pdp(
  model = mod,
  data = veteran,
  feature = "age",
  times = c(80, 160),
  method = "pdp+ice",
  grid.size = 8
)
plot_pdp(pdp_age, feature = "age", which = "per_time")

plot_pdp(pdp_age, feature = "age", which = "integrated", smooth = TRUE)
#> `geom_smooth()` using formula = 'y ~ x'

```
