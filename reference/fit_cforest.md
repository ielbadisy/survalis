# Fit a Conditional Inference Survival Forest

Fits a survival forest model using the
[`party::cforest()`](https://rdrr.io/pkg/party/man/cforest.html)
implementation of conditional inference trees for right-censored
survival data.

## Usage

``` r
fit_cforest(
  formula,
  data,
  teststat = "quad",
  testtype = "Univariate",
  mincriterion = 0,
  ntree = 500,
  mtry = 5,
  replace = TRUE,
  fraction = 0.632,
  ...
)
```

## Arguments

- formula:

  A [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) survival
  formula specifying the outcome and predictors.

- data:

  A data frame containing the variables in the model.

- teststat:

  Character string specifying the test statistic to use (default =
  `"quad"`). See
  [`party::cforest_control()`](https://rdrr.io/pkg/party/man/cforest_control.html).

- testtype:

  Character string specifying the type of test (default =
  `"Univariate"`).

- mincriterion:

  Numeric, the value of the test statistic that must be exceeded for a
  split to be performed (default = `0`).

- ntree:

  Integer, number of trees to grow (default = `500`).

- mtry:

  Integer, number of variables randomly selected at each node (default =
  `5`).

- replace:

  Logical, whether sampling of cases is with replacement (default =
  `TRUE`).

- fraction:

  Proportion of samples to draw if `replace = FALSE` (default =
  `0.632`).

- ...:

  Additional arguments passed to
  [`party::cforest_control()`](https://rdrr.io/pkg/party/man/cforest_control.html).

## Value

An object of class `"mlsurv_model"`, containing:

- model:

  The fitted `cforest` model.

- learner:

  Character string `"cforest"`.

- formula:

  The model formula.

- data:

  The training data.

## Details

This function wraps
[`party::cforest()`](https://rdrr.io/pkg/party/man/cforest.html) to fit
a conditional inference forest for survival analysis, returning an
`mlsurv_model` object compatible with
[`predict_cforest()`](https://elbadisyimad.com/survalis/reference/predict_cforest.md)
and the `survalis` framework.

## Examples

``` r
mod <- fit_cforest(Surv(time, status) ~ age + celltype + karno, data = veteran)
head(predict_cforest(mod, newdata = veteran[1:5, ], times = c(100, 200)))
#>       t=100     t=200
#> 1 0.5258457 0.2759512
#> 2 0.8285926 0.4212416
#> 3 0.5697253 0.4456722
#> 4 0.5800329 0.3240728
#> 5 0.8203375 0.3783192
```
