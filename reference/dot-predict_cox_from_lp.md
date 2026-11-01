# Reconstruct Cox-style survival curves from linear predictors

Internal helper shared by Cox-family learners to produce comparable
survival probabilities from a baseline cumulative hazard and
subject-specific linear predictors.

## Usage

``` r
.predict_cox_from_lp(lp_train, lp_new, y_train, times)
```

## Arguments

- lp_train:

  Numeric vector of linear predictors on the training data.

- lp_new:

  Numeric vector of linear predictors on `newdata`.

- y_train:

  A `Surv` object for the training data.

- times:

  Numeric evaluation times.

## Value

A standardized survival-probability data.frame.
