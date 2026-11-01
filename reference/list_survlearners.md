# List the survival learners

The ids accepted by the `model` argument of
[`fit()`](https://elbadisyimad.com/survalis/reference/fit.md),
[`evaluate()`](https://elbadisyimad.com/survalis/reference/evaluate.md),
[`compare()`](https://elbadisyimad.com/survalis/reference/compare.md),
[`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) and
[`estimate()`](https://elbadisyimad.com/survalis/reference/estimate.md).

## Usage

``` r
list_survlearners(has_tune = FALSE)
```

## Arguments

- has_tune:

  Logical (default `FALSE`). If `TRUE`, return only learners with a
  built-in tuning grid, so that
  [`tune()`](https://elbadisyimad.com/survalis/reference/tune.md) works
  without a `grid`.

## Value

A data frame with columns `learner` (the id) and `tunable` (whether a
built-in tuning grid exists).

## Examples

``` r
list_survlearners()
#>            learner tunable
#> 1            coxph   FALSE
#> 2            aalen   FALSE
#> 3           glmnet    TRUE
#> 4        selectcox    TRUE
#> 5           aftgee   FALSE
#> 6      flexsurvreg    TRUE
#> 7          bnnsurv    TRUE
#> 8            rpart    TRUE
#> 9             bart    TRUE
#> 10         xgboost    TRUE
#> 11        coxboost    TRUE
#> 12         fastgbm   FALSE
#> 13          ranger    TRUE
#> 14             rsf    TRUE
#> 15         cforest    TRUE
#> 16      blackboost    TRUE
#> 17         survsvm    TRUE
#> 18        densemlp   FALSE
#> 19            orsf    TRUE
#> 20          rpsurv   FALSE
#> 21          nbsurv   FALSE
#> 22 survmetalearner   FALSE
list_survlearners(has_tune = TRUE)
#>        learner tunable
#> 1       glmnet    TRUE
#> 2    selectcox    TRUE
#> 3  flexsurvreg    TRUE
#> 4      bnnsurv    TRUE
#> 5        rpart    TRUE
#> 6         bart    TRUE
#> 7      xgboost    TRUE
#> 8     coxboost    TRUE
#> 9       ranger    TRUE
#> 10         rsf    TRUE
#> 11     cforest    TRUE
#> 12  blackboost    TRUE
#> 13     survsvm    TRUE
#> 14        orsf    TRUE
```
