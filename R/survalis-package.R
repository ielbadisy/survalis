# R/survalis-package.R

#' SurvALIS: Interpretable Survival Machine Learning
#'
#' Core learners, tuning, evaluation, and interpretability utilities for survival analysis.
#'
#' @keywords internal
"_PACKAGE"



# Silence R CMD check NOTES for NSE variables used in dplyr/ggplot2 code
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(".", "effect", "observed_surv", "phi"))
}


#' Re-export Surv from survival
#' @importFrom survival Surv
#' @export
survival::Surv




#' @importFrom survival Surv survfit coxph aareg concordance
#' @importFrom flexsurv flexsurvreg
#' @importFrom randomForestSRC rfsrc
#' @importFrom partykit cforest
#' @importFrom timereg aalen cox.aalen
#' @importFrom mboost blackboost
#' @importFrom BART surv.bart
#' @importFrom xgboost xgb.train xgb.cv xgb.DMatrix
#' @importFrom gbm gbm
#' @importFrom nnet nnet
#' @importFrom rstpm2 stpm2 gsm nsx
#'
#' @importFrom stats model.frame model.matrix model.response delete.response
#'   terms formula update as.formula na.omit na.exclude na.pass predict
#'   runif rnorm rexp pnorm qnorm plogis sd quantile approx setNames
#'   coef predict.lm residuals lm glm binomial
#'
#' @importFrom utils tail head combn
#'
#' @importFrom methods is
#'
#' @importFrom dplyr select mutate arrange group_by summarise filter pull
#'   across everything all_of bind_cols bind_rows distinct
#'
#' @importFrom tidyr pivot_longer pivot_wider drop_na nest unnest
#'
#' @importFrom purrr map map2 pmap pmap_dbl imap_dbl map_dfr map_lgl
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_bar geom_histogram
#'   geom_boxplot labs theme_minimal theme element_text facet_wrap
#'   scale_color_manual scale_fill_manual
#'
#' @importFrom rlang .data enquo sym
#'
#' @importFrom future plan multisession sequential
#'
#' @importFrom parallel detectCores makeCluster stopCluster
#'
#' @importFrom progress progress_bar
#' @importFrom stats aggregate dist reorder reshape time var
#' @importFrom utils capture.output getFromNamespace
#' @importFrom ggplot2 geom_jitter geom_smooth geom_errorbar geom_abline coord_fixed coord_flip geom_col
#' @importFrom dplyr any_of pick
#' @importFrom stats pnorm
#' @importFrom pracma trapz
#' @importFrom data.table as.data.table setDT rbindlist
#' @importFrom gower gower_dist
#' @importFrom torch torch_tensor
#' @importFrom glmnet glmnet
#' @importFrom ggplot2 ggplot aes aes_string geom_line geom_point geom_col
#' @importFrom ggplot2 geom_errorbar geom_jitter geom_abline geom_smooth geom_vline
#' @importFrom ggplot2 labs theme_minimal coord_flip coord_fixed guides guide_legend
#' @importFrom ggplot2 position_dodge stat_summary scale_fill_manual scale_color_discrete
#' @importFrom ggplot2 geom_tile scale_fill_gradient
#' @importFrom glue glue


NULL




utils::globalVariables(c(
  # tidyselect / dplyr
  "metric","value","n","se","learner","scaled_importance","weight",
  "feature","feature1","feature2","feature_label","feature_value",
  "surv_prob","mean_pred_surv","lower_ci","upper_ci","calib_result",
  "bin","direction","count","time",
  # data.table-esque
  ".N",".id",
  # others seen in notes
  "splits","ale","integrated_ale","type"
))



