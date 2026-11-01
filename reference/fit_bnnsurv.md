# Fit a Bagged Nearest-Neighbour Survival Model

Fits a bagged ensemble of k-nearest-neighbour Kaplan-Meier learners
(Lowsky et al., 2013, bagged as in Wright's bnnSurvival). The engine is
implemented natively in survalis.

## Usage

``` r
fit_bnnsurv(
  formula,
  data,
  k = 5,
  num_base_learners = 10,
  num_features_per_base_learner = NULL,
  metric = c("mahalanobis", "euclidean"),
  weighting_function = function(x) x * 0 + 1,
  replace = TRUE,
  sample_fraction = NULL
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A data frame containing all variables referenced in `formula`.

- k:

  Integer, number of neighbours for each base learner. Default `5`.

- num_base_learners:

  Integer, number of base learners in the ensemble. Default `10`.

- num_features_per_base_learner:

  Integer or `NULL` (all predictors). Number of predictors drawn at
  random for each base learner.

- metric:

  Distance, `"mahalanobis"` (default) or `"euclidean"` (on standardized
  columns).

- weighting_function:

  Function mapping a vector of non-negative distances to non-negative
  neighbour weights. Default: constant weights.

- replace:

  Logical; sample rows with replacement for each base learner. Default
  `TRUE`.

- sample_fraction:

  Optional numeric in `(0, 1]`; fraction of rows drawn for each base
  learner. Defaults to `1` when `replace = TRUE` and `0.6321` otherwise.

## Value

An object of class `"mlsurv_model"` with elements:

- model:

  The fitted ensemble: design information, time grid and the base
  learners (row indices, columns and whitening maps).

- learner:

  Scalar `"bnnsurv"`.

- engine:

  Scalar `"survalis"`.

- formula, data:

  Inputs preserved for downstream use.

- time, status:

  Character names of the survival outcome fields.

## Method

Each base learner \\b = 1, \dots, B\\ draws a bootstrap (or subsample)
of the training rows and, optionally, a random subset of the predictors.
For a new observation \\x\\, the learner finds the \\k\\ training rows
closest to \\x\\ and returns their weighted Kaplan-Meier curve \$\$\hat
S_b(t \mid x) = \prod\_{t_j \le t} \left(1 - \frac{\sum_i w_i \delta_i
1\\T_i = t_j\\} {\sum_i w_i 1\\T_i \ge t_j\\}\right),\$\$ with weights
\\w_i = w(d(x, x_i))\\. The ensemble prediction is the average \\\hat
S(t \mid x) = B^{-1} \sum_b \hat S_b(t \mid x)\\, a non-increasing step
function in \\t\\. Training rows tied with the \\k\\-th nearest distance
share the remaining neighbour slots equally, so predictions do not
depend on the row order of the data.

The Mahalanobis distance uses the covariance \\\Sigma_b\\ of the base
learner's own sample: \\d(x, y)^2 = (x - y)^\top \Sigma_b^{+} (x - y)\\,
where \\\Sigma_b^{+}\\ is the Moore-Penrose pseudo-inverse. Constant and
collinear columns (for example a dummy that never varies in a bootstrap
sample) therefore contribute nothing instead of breaking the inversion.
`metric = "euclidean"` uses the Euclidean distance on standardized
columns.

Differences from bnnSurvival: factors are one-hot encoded (not coded as
integer labels) and new data reuse the training levels; feature subsets
are drawn from all predictors (a factor counts as one predictor);
constant columns are detected from the variances; prediction runs in the
current process (no forking); and predictions at arbitrary times are
read off the step function rather than linearly interpolated. With
numeric predictors, all predictors per learner, constant weights and no
tied distances, the two engines give the same curves for the same seed.

## References

Lowsky DJ, Ding Y, Lee DKK, McCulloch CE, Ross LF, Thistlethwaite JR,
Zenios SA (2013). A K-nearest neighbors survival probability prediction
method. *Statistics in Medicine*, 32(12), 2062-2069.
[doi:10.1002/sim.5673](https://doi.org/10.1002/sim.5673)

## See also

[`predict_bnnsurv()`](https://elbadisyimad.com/survalis/reference/predict_bnnsurv.md),
[`tune_bnnsurv()`](https://elbadisyimad.com/survalis/reference/tune_bnnsurv.md)

## Examples

``` r
mod <- fit_bnnsurv(
  Surv(time, status) ~ age + karno + diagtime + prior,
  data = veteran
)
head(predict_bnnsurv(mod, newdata = veteran[1:5, ], times = c(50, 100)))
#>   t=50 t=100
#> 1 0.82  0.28
#> 2 0.98  0.62
#> 3 0.72  0.34
#> 4 0.78  0.78
#> 5 0.94  0.94
```
