# Standardized RMST Estimate or Causal-Style RMST Contrast via G-Computation

Estimates, from a fitted survalis model (`fit_*()` output), either (i)
the standardized marginal restricted mean survival time (RMST) at a
fixed horizon, or (ii) a covariate-adjusted contrast in RMST between the
two levels of a binary grouping variable via g-computation – with
percentile bootstrap inference adapted from the tvrmst package's
dynamic-RMST bootstrap design. The two variants are dispatched by the
`trt_col` argument (`NULL` = marginal estimate, matching the
`group = NULL` dispatch already used by
[`plot_survmat`](https://elbadisyimad.com/survalis/reference/plot_survmat.md)),
rather than by separate function names.

## Usage

``` r
estimate_rmst(
  model,
  data,
  tau,
  times,
  trt_col = NULL,
  R = 200,
  conf = 0.95,
  seed = NULL
)
```

## Arguments

- model:

  An `mlsurv_model` object from a `fit_*()` learner (must expose
  `model$learner` so `predict_<learner>()` can be dispatched, as in
  [`compute_calibration`](https://elbadisyimad.com/survalis/reference/compute_calibration.md)).

- data:

  A data.frame of covariates to standardize over.

- tau:

  Positive numeric restriction time for RMST.

- times:

  Numeric vector of evaluation times passed to `predict_<learner>()`;
  must be dense enough to integrate up to `tau` (see
  [`survmat_to_rmst`](https://elbadisyimad.com/survalis/reference/survmat_to_rmst.md)).

- trt_col:

  Optional character; name of a binary grouping/exposure column in
  `data`. If `NULL` (default), `estimate_rmst()` returns the
  standardized marginal RMST at `tau` for `data` as given. If supplied,
  its two observed values are used as counterfactual levels (the larger
  value is the "treated"/A=1 arm) and the function returns the
  standardized contrast between them.

- R:

  Number of bootstrap replicates. Default `200`.

- conf:

  Confidence level in (0, 1). Default `0.95`.

- seed:

  Optional RNG seed for reproducibility.

## Value

A list of class `"rmst_estimate"` with elements `estimate`, `lo`, `hi`,
`conf`, `tau`, `trt_col` (`NULL` in marginal mode), `levels` (the two
levels compared in contrast mode, or `NULL`), `mode` (`"marginal"` or
`"contrast"`), and `boot` (the vector of bootstrap replicate estimates).

## Details

**Marginal mode (`trt_col = NULL`).** Individual RMSTs
\\\mathrm{RMST}\_i(\tau)\\ are obtained via
[`survmat_to_rmst`](https://elbadisyimad.com/survalis/reference/survmat_to_rmst.md)
from \\\hat S(t \mid X_i)\\, and the standardized estimate is
\\\widehat{\mathrm{RMST}}(\tau) = n^{-1}\sum_i \mathrm{RMST}\_i(\tau)\\.

**Contrast mode (`trt_col` supplied).** Counterfactual survival curves
\\\hat S(t \mid X_i, A{=}1)\\ and \\\hat S(t \mid X_i, A{=}0)\\ are
obtained by setting `trt_col` to each level for the entire data set and
predicting from the fixed, already-fitted `model`
(g-computation/standardization). The standardized contrast is \$\$
\widehat{\Delta}\_{\mathrm{RMST}}(\tau) = \frac{1}{n}\sum\_{i=1}^n
\mathrm{RMST}\_i(\tau, A{=}1) - \frac{1}{n}\sum\_{i=1}^n
\mathrm{RMST}\_i(\tau, A{=}0). \$\$

In both modes, uncertainty is quantified by a nonparametric percentile
bootstrap over subjects (rows of `data`), holding `model` fixed (i.e.
conditional on the fitted ensemble, as in tvrmst's `boot_rmst_delta()`):
each replicate resamples rows with replacement and recomputes the
estimate.

When `trt_col` corresponds to a manipulable intervention and standard
causal identification assumptions hold,
\\\widehat{\Delta}\_{\mathrm{RMST}}(\tau)\\ admits a causal
interpretation; for non-manipulable attributes (e.g. sex, biomarker
status) it should be read as a covariate-adjusted marginal contrast
rather than a causal effect.

## References

Royston P, Parmar MKB (2013). "Restricted Mean Survival Time: An
Alternative to the Hazard Ratio for the Design and Analysis of
Randomized Trials with a Time-to-Event Outcome." *BMC Medical Research
Methodology*, 13, 152.

## See also

[`survmat_to_rmst()`](https://elbadisyimad.com/survalis/reference/survmat_to_rmst.md),
[`plot_rmst()`](https://elbadisyimad.com/survalis/reference/plot_rmst.md)

## Examples

``` r
# \donttest{
mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
times <- seq(10, 300, by = 10)

marg <- estimate_rmst(mod, data = veteran, tau = 200, times = times, R = 20, seed = 1)
marg$estimate
#> [1] 89.52165

contrast <- estimate_rmst(mod, data = veteran, tau = 200, times = times,
                          trt_col = "trt", R = 20, seed = 1)
contrast$estimate
#> [1] -9.565834
# }
```
