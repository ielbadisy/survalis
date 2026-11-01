# Local Surrogate Explanation for Survival Predictions (LIME-style)

Builds a sparse, locally weighted linear surrogate model to explain a
fitted survival model's prediction at a specific target time for a
single instance. Categorical features are binarized relative to the
instance of interest; local weights are computed from feature-space
proximity (Gower by default); and a penalized (lasso) or unpenalized
linear model is fit to approximate the model's predicted survival
probability at `target_time`.

## Usage

``` r
compute_surrogate(
  model,
  newdata,
  baseline_data,
  times,
  target_time,
  k = 5,
  dist.fun = "gower",
  gower.power = 5,
  kernel.width = NULL,
  penalized = TRUE,
  exclude = NULL
)
```

## Arguments

- model:

  An `mlsurv_model` fitted via `fit_*()`, containing a valid `$learner`
  so that `predict_<learner>` can be dispatched.

- newdata:

  A one-row data frame: the instance to explain (must have the same
  predictors as `baseline_data`).

- baseline_data:

  A data frame used to define the local neighborhood and to fit the
  surrogate (typically the model's training data).

- times:

  Numeric vector of times passed to the prediction function; must
  include `target_time` (nearest value is used).

- target_time:

  Numeric time at which to explain the prediction.

- k:

  Desired number of non-zero coefficients for the penalized surrogate
  (used only when `penalized = TRUE`); acts as a target sparsity.

- dist.fun:

  Distance function for locality weighting. Default `"gower"` (via
  gower); otherwise any method accepted by
  [`dist`](https://rdrr.io/r/stats/dist.html) (requires `kernel.width`).

- gower.power:

  Power applied to \\(1 - \text{gower\\dist})\\ when
  `dist.fun = "gower"` to control locality sharpness (default `5`).

- kernel.width:

  Positive numeric bandwidth used for non-Gower kernels (Gaussian
  weighting on pairwise distances).

- penalized:

  Logical; if `TRUE` (default) fits a lasso surrogate via glmnet;
  otherwise fits a weighted least squares
  [`lm()`](https://rdrr.io/r/stats/lm.html).

- exclude:

  Optional character vector of column names to exclude from the
  surrogate (in addition to survival outcome columns).

## Value

A data frame with one row per selected feature containing:

- feature:

  Feature name (after recoding).

- feature_value:

  Value of the feature in `newdata`.

- effect:

  Local contribution \\\beta_j \cdot x_j\\ at `target_time`.

- target_time:

  The explanation time (copied for plotting).

Rows are ordered by decreasing \\\|effect\|\\.

## Details

**Target:** The surrogate approximates \\S(t^{\*} \mid x)\\ returned by
the underlying model at `target_time`. The nearest column in `times` is
used.

**Weights:** Locality weights are computed from distances between
`baseline_data` rows and `newdata`. With `"gower"`, weights are \\(1 -
\mathrm{gower\\dist})^{\mathrm{gower.power}}\\; otherwise Gaussian
\\\exp\\-d^2 / \mathrm{kernel.width}^2\\^{1/2}\\.

**Sparsity:** When `penalized=TRUE`, a glmnet path is fit and the
solution with degrees of freedom closest to `k` is selected (preferring
exactly `k` if available).

## Examples

``` r
# \donttest{
mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
local_expl <- compute_surrogate(
  model = mod,
  newdata = veteran[2, , drop = FALSE],
  baseline_data = veteran,
  times = c(80, 160),
  target_time = 80,
  k = 3
)
head(local_expl)
#>    feature feature_value     effect target_time
#> 1    karno            70 0.65901926          80
#> 2 celltype      squamous 0.18321600          80
#> 3      age            64 0.04473744          80
# }
```
