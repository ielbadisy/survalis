
learners <- tibble::tibble(
  learner_name = c("aareg", "aftgee", "mboost", "ranger", "orsf", "glmnet", "flexsurv", "bnnsurv", "coxaalen"),
  fit_fun = list(fit_aareg, fit_aftgee, fit_mboost, fit_ranger, fit_orsf, fit_glmnet, fit_flexsurv, fit_bnnsurv, fit_coxaalen),
  pred_fun = list(predict_aareg, predict_aftgee, predict_mboost, predict_ranger, predict_orsf, predict_glmnet, predict_flexsurv, predict_bnnsurv, predict_coxaalen)
)


results_by_learner <- purrr::pmap_dfr(
  learners,
  function(learner_name, fit_fun, pred_fun) {
    message("Testing: ", learner_name)
    tryCatch({
      res <- cv_survlearner(
        formula = Surv(Time, status) ~ x1 + x2,
        data = dat,
        fit_fun = fit_fun,
        pred_fun = pred_fun,
        times = c(10, 30, 50),
        metrics = c("cindex", "ibs"),
        folds = 6,
        seed = 123
      )
      dplyr::mutate(res, learner_name = learner_name)
    }, error = function(e) {
      warning("Learner '", learner_name, "' failed: ", conditionMessage(e))
      NULL
    })
  }
)
