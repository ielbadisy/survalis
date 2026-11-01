# List tunable survival learners

Filters learners to only those that provide a `tune_*` function in
addition to `fit_*` and `predict_*`.

## Usage

``` r
list_tunable_survlearners()
```

## Value

A base `data.frame` like
[`list_survlearners()`](https://elbadisyimad.com/survalis/reference/list_survlearners.md)
but containing only rows where `tune` is not `NA`.

## Examples

``` r
list_tunable_survlearners()
#>        learner             fit             predict             tune has_fit
#> 3       glmnet      fit_glmnet      predict_glmnet      tune_glmnet    TRUE
#> 4    selectcox   fit_selectcox   predict_selectcox   tune_selectcox    TRUE
#> 6  flexsurvreg fit_flexsurvreg predict_flexsurvreg tune_flexsurvreg    TRUE
#> 7      bnnsurv     fit_bnnsurv     predict_bnnsurv     tune_bnnsurv    TRUE
#> 8        rpart       fit_rpart       predict_rpart       tune_rpart    TRUE
#> 9         bart        fit_bart        predict_bart        tune_bart    TRUE
#> 10     xgboost     fit_xgboost     predict_xgboost     tune_xgboost    TRUE
#> 11    coxboost    fit_coxboost    predict_coxboost    tune_coxboost    TRUE
#> 13      ranger      fit_ranger      predict_ranger      tune_ranger    TRUE
#> 14         rsf         fit_rsf         predict_rsf         tune_rsf    TRUE
#> 15     cforest     fit_cforest     predict_cforest     tune_cforest    TRUE
#> 16  blackboost  fit_blackboost  predict_blackboost  tune_blackboost    TRUE
#> 17     survsvm     fit_survsvm     predict_survsvm     tune_survsvm    TRUE
#> 19        orsf        fit_orsf        predict_orsf        tune_orsf    TRUE
#>    has_predict has_tune available
#> 3         TRUE     TRUE      TRUE
#> 4         TRUE     TRUE      TRUE
#> 6         TRUE     TRUE      TRUE
#> 7         TRUE     TRUE      TRUE
#> 8         TRUE     TRUE      TRUE
#> 9         TRUE     TRUE      TRUE
#> 10        TRUE     TRUE      TRUE
#> 11        TRUE     TRUE      TRUE
#> 13        TRUE     TRUE      TRUE
#> 14        TRUE     TRUE      TRUE
#> 15        TRUE     TRUE      TRUE
#> 16        TRUE     TRUE      TRUE
#> 17        TRUE     TRUE      TRUE
#> 19        TRUE     TRUE      TRUE
```
