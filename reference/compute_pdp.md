# Partial Dependence and ICE for Survival Predictions

Computes partial dependence (PDP) and/or individual conditional
expectation (ICE) curves of predicted survival probabilities for a
single feature at one or more evaluation times. Works with any learner
fitted via `fit_*()` that exposes a matching `predict_*()` method
returning survival probabilities.

## Usage

``` r
compute_pdp(model, data, feature, times, method = "pdp+ice", grid.size = 20)
```

## Arguments

- model:

  An `mlsurv_model` created by a `fit_*()` function; must contain a
  valid `learner` so the appropriate `predict_<learner>()` can be
  dispatched.

- data:

  A data frame used to construct PDP/ICE profiles (typically the
  training data).

- feature:

  Character scalar; the feature name to analyze (numeric or
  categorical).

- times:

  Numeric vector of evaluation times at which survival probabilities are
  computed.

- method:

  One of `"pdp"`, `"ice"`, or `"pdp+ice"` (default). Controls which
  profiles are produced.

- grid.size:

  Integer number of grid points for numeric features (default 20).
  Ignored for categorical features (levels are used).

## Value

A list with elements:

- results:

  Long data frame with columns: `surv_prob`, `time`, `type` (pdp/ice),
  `.id` (row id for ICE), and the analyzed `feature`.

- pdp_integrated:

  (Optional) Data frame with `feature` and `integrated_surv`
  (time-integrated PDP), present only if `length(times) > 1` and PDP was
  requested.

## Details

For numeric features, a regular grid over the observed range is used;
for categorical features, all observed levels are used. For each grid
value, predictions are made for every row of `data` with the feature
forced to the grid value, yielding ICE curves per row and PDP as the
average across rows.

If multiple `times` are supplied, outputs are stacked in long format
with a `time` column; an additional integrated PDP is computed via
trapezoidal rule when more than one time is provided.

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
head(pdp_age$results)
#>            age .id time surv_prob type
#> ice.1 34.00000   1   80 0.4962103  ice
#> ice.2 40.71429   1   80 0.5051970  ice
#> ice.3 47.42857   1   80 0.5141100  ice
#> ice.4 54.14286   1   80 0.5229459  ice
#> ice.5 60.85714   1   80 0.5317016  ice
#> ice.6 67.57143   1   80 0.5403740  ice
```
