# Fit a Parametric Survival Regression Model Using flexsurvreg

This function fits a fully parametric survival regression model using
the flexsurv package. It supports a variety of parametric distributions
(e.g., Weibull, exponential, log-normal) and returns a model object
compatible with the `mlsurv_model` interface for downstream survival
predictions.

## Usage

``` r
fit_flexsurvreg(formula, data, dist = "weibull", ...)
```

## Arguments

- formula:

  A [`Surv`](https://rdrr.io/pkg/survival/man/Surv.html)-based survival
  formula of the form `Surv(time, status) ~ predictors`.

- data:

  A data frame containing the variables in the model.

- dist:

  Character string specifying the parametric distribution to use
  (default = `"weibull"`). See
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)
  for available options.

- ...:

  Additional arguments passed to
  [`flexsurvreg`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md).

## Value

An object of class `mlsurv_model` containing:

- `model` - the fitted `flexsurvreg` model

- `learner` - string identifier `"flexsurvreg"`

- `formula` - the model formula

- `data` - the training data

- `time` - the name of the time-to-event column

- `status` - the name of the event indicator column

## Examples

``` r

mod_flex <- fit_flexsurvreg(Surv(time, status) ~ age + celltype + karno,
                            data = veteran,
                            dist = "weibull")
```
