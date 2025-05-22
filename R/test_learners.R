learners <- tibble::tibble(
  learner_name = c("aareg", "aftgee", "bart", "mboost"),
  fit_fun = list(fit_aareg, fit_aftgee, fit_bart, fit_mboost),
  pred_fun = list(predict_aareg, predict_aftgee, predict_bart, predict_mboost)
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
        times = c(10, 20, 40),
        metrics = c("cindex", "ibs"),
        folds = 5,
        seed = 123
      )
      dplyr::mutate(res, learner_name = learner_name)
    }, error = function(e) {
      warning("Learner '", learner_name, "' failed: ", conditionMessage(e))
      NULL
    })
  }
)
