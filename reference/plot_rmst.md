# RMST Estimate or Contrast Curve Over a Grid of Restriction Times

Evaluates
[`estimate_rmst`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md)
across a grid of restriction times \\\tau\\ and plots the resulting
curve with percentile bootstrap confidence bands. As in
[`estimate_rmst`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md),
marginal vs. contrast mode is dispatched by `trt_col` rather than by
separate function names.

## Usage

``` r
plot_rmst(
  model,
  data,
  tau_seq,
  times,
  trt_col = NULL,
  R = 100,
  conf = 0.95,
  seed = NULL,
  title
)
```

## Arguments

- model:

  An `mlsurv_model` object from a `fit_*()` learner.

- data:

  A data.frame of covariates.

- tau_seq:

  Numeric vector of restriction times to evaluate.

- times:

  Numeric vector of evaluation times passed to `predict_<learner>()`.

- trt_col:

  Optional character; name of a binary grouping/exposure column. `NULL`
  (default) plots the marginal RMST(tau) curve; supplying it plots the
  RMST contrast curve.

- R:

  Number of bootstrap replicates per `tau`. Default `100`.

- conf:

  Confidence level in (0, 1). Default `0.95`.

- seed:

  Optional RNG seed.

- title:

  Plot title. If missing, an automatically generated title is used. Pass
  `NULL` to omit the title entirely.

## Value

A ggplot2 object with an attribute `"curve"` holding the underlying
`data.frame` (`tau`, `estimate`, `lo`, `hi`).

## See also

[`estimate_rmst()`](https://elbadisyimad.com/survalis/reference/estimate_rmst.md)

## Examples

``` r
# \donttest{
mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
times <- seq(10, 300, by = 10)
plot_rmst(mod, data = veteran, tau_seq = seq(50, 250, by = 50),
                times = times, trt_col = "trt", R = 10, seed = 1)

# }
```
