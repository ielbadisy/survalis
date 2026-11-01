# Standardised causal-style contrast for a binary treatment

Estimates a covariate-adjusted marginal contrast between the two levels
of a binary `treatment` by G-computation (outcome standardisation): a
survival model is fitted for
`Surv(time, status) ~ treatment * (confounders)`, every subject's
outcome is predicted under both treatment levels, and the counterfactual
mean survival curves are contrasted.

## Usage

``` r
estimate(
  formula,
  data,
  model,
  estimand = c("RMST_diff", "RMST_ratio", "survival_diff", "survival_ratio"),
  treatment,
  tau,
  times = NULL,
  spec = list(),
  R = 200,
  conf = 0.95,
  refit = TRUE,
  seed = NULL,
  ...
)
```

## Arguments

- formula:

  A survival formula whose right-hand side lists the **confounders**
  (the treatment is added by the function, with treatment-by-confounder
  interactions).

- data:

  A data frame.

- model:

  A learner id (`list_survlearners()$learner`) or a custom learner (see
  [`fit()`](https://elbadisyimad.com/survalis/reference/fit.md)).

- estimand:

  One of `"RMST_diff"` (default), `"RMST_ratio"`, `"survival_diff"`,
  `"survival_ratio"`.

- treatment:

  Name of a binary column in `data`. The larger value is the "treated"
  (A = 1) level.

- tau:

  Horizon: the RMST restriction time for the RMST estimands, or the
  evaluation time for the survival estimands.

- times:

  Prediction grid; `NULL` picks up to 50 event-time quantiles (5th to
  95th percentile). Must reach `tau`.

- spec:

  Engine arguments for the learner (named list).

- R:

  Bootstrap replicates (default `200`).

- conf:

  Confidence level (default `0.95`).

- refit:

  Logical; refit the model inside each bootstrap replicate (default
  `TRUE`).

- seed:

  Optional RNG seed.

- ...:

  Unused.

## Value

A `survalis_estimate` object with `estimand`, `estimate`, `conf.low`,
`conf.high`, `conf`, `treatment`, `levels`, `tau`, `R`, `refit`, `boot`
(replicate estimates) and `curves` (counterfactual mean survival for
each arm on `times`).

## Details

Uncertainty is a nonparametric bootstrap over subjects; by default the
model is refitted in each replicate (`refit = TRUE`), which propagates
model estimation error. Set `refit = FALSE` for a faster interval that
is conditional on the fitted model.

## Interpretation

The estimate has a causal reading only under standard identification
assumptions: consistency, no unmeasured confounding given the `formula`
covariates, positivity, no interference, and a correctly specified
outcome model. For a non-manipulable attribute it is a
covariate-adjusted marginal contrast, not a causal effect. This is a
plug-in G-formula estimator; it is not doubly robust.

## Examples

``` r
est <- estimate(Surv(time, status) ~ age + karno + celltype, veteran,
                model = "coxph", estimand = "RMST_diff",
                treatment = "trt", tau = 300, R = 50, seed = 1)
est
#> <survalis_estimate> RMST difference at tau = 300  (model: coxph)
#>   treatment: trt  (2 vs 1)
#>   estimate: -14.3905   95% CI [-40.6135, 7.8719]   (50 bootstrap, refit)
#>   plug-in G-formula; causal reading needs exchangeability, positivity,
#>   consistency, no interference and a correct outcome model.
```
