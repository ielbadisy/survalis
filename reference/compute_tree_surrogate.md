# Compute Tree-Based Surrogate Model for Survival Predictions

Fits decision tree surrogate models to approximate the predictions of a
fitted survival model at one or more evaluation times. This allows users
to gain interpretable, rule-based approximations of complex survival
models.

## Usage

``` r
compute_tree_surrogate(model, data, times, minsplit = 10, cp = 0.01)
```

## Arguments

- model:

  A fitted survival model object created with a `fit_*()` function. Must
  include a valid `learner` field corresponding to a `predict_*()`
  method.

- data:

  A data frame containing predictor variables (and optional survival
  outcome columns).

- times:

  A numeric vector of evaluation times at which to approximate model
  predictions. Must contain at least one value.

- minsplit:

  Minimum number of observations required to attempt a split in the
  surrogate tree. Passed to
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html).

- cp:

  Complexity parameter for the surrogate tree. Passed to
  [`rpart::rpart.control()`](https://rdrr.io/pkg/rpart/man/rpart.control.html).

## Value

An object of class `"tree_surrogate"`, containing:

- `times`: the evaluation times.

- `results`: a list with one element per time, each containing:

  - `tree`: the fitted `rpart` object.

  - `r_squared`: the R of the surrogate model vs. the original
    predictions.

  - `split_count`: a table of feature split counts.

- `dynamic`: logical indicating if more than one time was used.

## Details

For each evaluation time, the function:

1.  Predicts survival probabilities from the fitted model.

2.  Excludes survival outcome columns (`time`, `status`, `event`) from
    the predictors.

3.  Fits a decision tree to approximate the predicted probabilities.

4.  Computes the R between the model predictions and the surrogate
    predictions.

5.  Counts the number of splits per feature.

If multiple `times` are provided, results are stored for each time
point.

## Examples

``` r
if (FALSE) { # \dontrun{
mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
tree_ranger <- compute_tree_surrogate(
  model = mod_ranger,
  data = veteran,
  times = c(100, 200, 300)
)
} # }
```
