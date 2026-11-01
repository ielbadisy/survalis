# Compute Feature Interactions for Survival Predictions

Estimates global and time-varying interaction strengths of model
predictions, using a Friedman-H style decomposition adapted to survival
partial dependence. Works with any `mlsurv_model` that has a matching
`predict_*()` method returning survival probabilities.

## Usage

``` r
compute_interactions(
  model,
  data,
  times,
  target_time = NULL,
  features = NULL,
  type = c("1way", "heatmap", "time"),
  grid.size = 30,
  batch.size = 100
)
```

## Arguments

- model:

  An `mlsurv_model` created via a `fit_*()` function; must include a
  valid `learner` so the appropriate `predict_<learner>()` can be
  dispatched.

- data:

  Data frame used to probe the model (typically training data).

- times:

  Numeric vector of evaluation times used for prediction.

- target_time:

  Single time at which to quantify interactions for `type = "1way"` or
  `type = "heatmap"`. Required for these types.

- features:

  Optional character vector of feature names to evaluate. Defaults to
  all predictors in `data` except the outcome columns.

- type:

  One of:

  - `"1way"` – per-feature Friedman-H interaction strength vs. all
    others at `target_time`;

  - `"heatmap"` – pairwise feature interaction matrix at `target_time`;

  - `"time"` – per-feature interaction strength across `times`.

- grid.size:

  Integer; number of random grid values / replicates used for Monte
  Carlo marginalization (default 30).

- batch.size:

  Reserved for future batching support (currently unused).

## Value

A data frame whose structure depends on `type`:

- `"1way"`: columns `feature`, `interaction`.

- `"heatmap"`: columns `feature1`, `feature2`, `interaction` (symmetric
  with zeros on the diagonal).

- `"time"`: columns `feature`, `time`, `interaction`.

## Details

For a target time \\t^{\*}\\, let \\f(x)\\ be the predicted survival
probability. For feature \\j\\, we approximate a decomposition: \$\$
f(x) \approx f_j(x_j) + f\_{-j}(x\_{-j}) \$\$ by Monte Carlo
marginalization over subsets of features using random sampling from the
empirical distribution in `data`. The reported interaction strength is:
\$\$ H_j(t^{\*}) = \sqrt{ \frac{\sum ( f(x) - \\ \tilde f_j(x_j) +
\tilde f\_{-j}(x\_{-j}) \\)^2 }{ \sum f(x)^2 } } , \$\$ clipped to \[0,
1\]. Pairwise heatmaps are computed analogously for \\(j, k)\\.

Larger values indicate stronger non-additivity (interaction) involving
the feature(s). The `"time"` mode repeats the computation across all
`times` to show dynamics.

## References

Friedman, J. H., and Popescu, B. E. (2008). Predictive learning via rule
ensembles. *Annals of Applied Statistics*. (Friedman's \\H\\ interaction
measure.)

## Examples

``` r
mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
times <- c(80, 160)
compute_interactions(
  model = mod,
  data = veteran,
  times = times,
  target_time = 80,
  features = c("age", "karno"),
  type = "1way",
  grid.size = 6
)
#>   feature interaction
#> 1     age   0.8789678
#> 2   karno   1.0000000
```
