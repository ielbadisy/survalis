# Fit a Survival Tree Model using `rpart`

Fits a survival tree using the rpart package with an exponential
splitting rule. This learner is compatible with the `mlsurv_model`
interface and can be used with the
[`cv_survlearner()`](https://elbadisyimad.com/survalis/reference/cv_survlearner.md)
and `tune_*()` functions in the survalis framework.

## Usage

``` r
fit_rpart(
  formula,
  data,
  minsplit = 20,
  minbucket = round(minsplit/3),
  cp = 0.01,
  maxcompete = 4,
  maxsurrogate = 5,
  usesurrogate = 2,
  xval = 10,
  surrogatestyle = 0,
  maxdepth = 30
)
```

## Arguments

- formula:

  A survival formula of the form `Surv(time, status) ~ predictors`.

- data:

  A `data.frame` containing the variables in the model.

- minsplit:

  Minimum number of observations that must exist in a node in order for
  a split to be attempted. Default is `20`.

- minbucket:

  Minimum number of observations in any terminal node. Default is
  `round(minsplit / 3)`.

- cp:

  Complexity parameter for pruning. Default is `0.01`.

- maxcompete:

  Number of competitor splits retained in the output. Default is `4`.

- maxsurrogate:

  Number of surrogate splits retained in the output. Default is `5`.

- usesurrogate:

  How surrogates are used in the splitting process. Default is `2`.

- xval:

  Number of cross-validations to perform in `rpart`. Default is `10`.

- surrogatestyle:

  Controls selection of surrogate splits. Default is `0`.

- maxdepth:

  Maximum depth of any node of the final tree. Default is `30`.

## Value

An object of class `"mlsurv_model"` containing:

- `model` – fitted `rpart` survival tree object

- `learner` – `"rpart"`

- `formula`, `data`, `time`, and `status`

## Details

This function fixes the `method = "exp"` argument to fit an exponential
splitting survival tree, which models the hazard function assuming
exponential survival within terminal nodes.

## References

Therneau TM, Atkinson EJ. (2019). *An Introduction to Recursive
Partitioning Using the RPART Routines*. Mayo Clinic.

## Examples

``` r
mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred_rpart <- predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300))
head(pred_rpart)
#>       t=100     t=200      t=300
#> 1 0.4062035 0.1650013 0.06702412
#> 2 0.6080413 0.3697142 0.22480154
#> 3 0.4062035 0.1650013 0.06702412
#> 4 0.4062035 0.1650013 0.06702412
#> 5 0.6080413 0.3697142 0.22480154
```
