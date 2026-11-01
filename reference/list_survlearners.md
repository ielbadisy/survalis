# List survival learners available in survalis

Returns a table of known survival learners, showing the expected
`fit_*`, `predict_*`, and (if present) `tune_*` functions, plus booleans
indicating which functions are available.

## Usage

``` r
list_survlearners(has_tune = FALSE)
```

## Arguments

- has_tune:

  Logical (default `FALSE`). If `TRUE`, return only learners that have a
  corresponding `tune_*` function. If `FALSE`, return all.

## Value

A data.frame with columns:

- `learner` – learner id (e.g., "ranger")

- `fit`, `predict`, `tune` – function names (tune may be `NA`)

- `has_fit`, `has_predict`, `has_tune` – logical flags

- `available` – `has_fit & has_predict`

## Examples

``` r
list_survlearners()                   # all learners (default)
#>            learner                 fit                 predict             tune
#> 1            coxph           fit_coxph           predict_coxph             <NA>
#> 2            aalen           fit_aalen           predict_aalen             <NA>
#> 3           glmnet          fit_glmnet          predict_glmnet      tune_glmnet
#> 4        selectcox       fit_selectcox       predict_selectcox   tune_selectcox
#> 5           aftgee          fit_aftgee          predict_aftgee             <NA>
#> 6      flexsurvreg     fit_flexsurvreg     predict_flexsurvreg tune_flexsurvreg
#> 7          bnnsurv         fit_bnnsurv         predict_bnnsurv     tune_bnnsurv
#> 8            rpart           fit_rpart           predict_rpart       tune_rpart
#> 9             bart            fit_bart            predict_bart        tune_bart
#> 10         xgboost         fit_xgboost         predict_xgboost     tune_xgboost
#> 11        coxboost        fit_coxboost        predict_coxboost    tune_coxboost
#> 12         fastgbm         fit_fastgbm         predict_fastgbm             <NA>
#> 13          ranger          fit_ranger          predict_ranger      tune_ranger
#> 14             rsf             fit_rsf             predict_rsf         tune_rsf
#> 15         cforest         fit_cforest         predict_cforest     tune_cforest
#> 16      blackboost      fit_blackboost      predict_blackboost  tune_blackboost
#> 17         survsvm         fit_survsvm         predict_survsvm     tune_survsvm
#> 18        densemlp        fit_densemlp        predict_densemlp             <NA>
#> 19            orsf            fit_orsf            predict_orsf        tune_orsf
#> 20          rpsurv          fit_rpsurv          predict_rpsurv             <NA>
#> 21          nbsurv          fit_nbsurv          predict_nbsurv             <NA>
#> 22 survmetalearner fit_survmetalearner predict_survmetalearner             <NA>
#>    has_fit has_predict has_tune available
#> 1     TRUE        TRUE    FALSE      TRUE
#> 2     TRUE        TRUE    FALSE      TRUE
#> 3     TRUE        TRUE     TRUE      TRUE
#> 4     TRUE        TRUE     TRUE      TRUE
#> 5     TRUE        TRUE    FALSE      TRUE
#> 6     TRUE        TRUE     TRUE      TRUE
#> 7     TRUE        TRUE     TRUE      TRUE
#> 8     TRUE        TRUE     TRUE      TRUE
#> 9     TRUE        TRUE     TRUE      TRUE
#> 10    TRUE        TRUE     TRUE      TRUE
#> 11    TRUE        TRUE     TRUE      TRUE
#> 12    TRUE        TRUE    FALSE      TRUE
#> 13    TRUE        TRUE     TRUE      TRUE
#> 14    TRUE        TRUE     TRUE      TRUE
#> 15    TRUE        TRUE     TRUE      TRUE
#> 16    TRUE        TRUE     TRUE      TRUE
#> 17    TRUE        TRUE     TRUE      TRUE
#> 18    TRUE        TRUE    FALSE      TRUE
#> 19    TRUE        TRUE     TRUE      TRUE
#> 20    TRUE        TRUE    FALSE      TRUE
#> 21    TRUE        TRUE    FALSE      TRUE
#> 22    TRUE        TRUE    FALSE      TRUE
list_survlearners(has_tune = TRUE)    # only tunable learners
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
