fit_selectcox <- function(formula, data, rule = "aic") {
  stopifnot(requireNamespace("pec", quietly = TRUE))

  model <- pec::selectCox(formula = formula, data = data, rule = rule)

  structure(list(
    model = model,
    learner = "selectcox",
    formula = formula,
    data = data,
    rule = rule
  ), class = "mlsurv_model", engine = "selectCox")
}

predict_selectcox <- function(object, newdata, times) {
  stopifnot(requireNamespace("pec", quietly = TRUE))

  pred <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  pred <- as.matrix(pred)

  colnames(pred) <- paste0("t=", times)
  rownames(pred) <- paste0("ID_", seq_len(nrow(pred)))

  as.data.frame(pred)
}



tune_selectcox <- function(formula, data, times,
                           rules = c("aic", "p"),
                           metrics = c("cindex", "ibs"),
                           folds = 5, seed = 123, ...) {

  purrr::map_dfr(rules, function(rule) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_selectcox,
      pred_fun = predict_selectcox,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      rule = rule,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(rule = rule) |>
      dplyr::bind_cols(tidyr::pivot_wider(summary[, c("metric", "mean")], names_from = metric, values_from = mean))
  }) |> arrange(dplyr::across(any_of(metrics[1])))
}

res_selectcox <- tune_selectcox(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  rules = c("aic", "p"),
  metrics = c("cindex", "ibs", "ise"),
  folds = 3
)

print(res_selectcox)

