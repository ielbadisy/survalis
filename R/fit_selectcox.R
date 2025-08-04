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


  if (!is.null(object$learner) && object$learner != "selectcox") {
    warning("Object passed to predict_selectcox() may not come from fit_selectcox().")
  }


  stopifnot(requireNamespace("pec", quietly = TRUE))

  pred <- pec::predictSurvProb(object$model, newdata = newdata, times = times)
  survmat <- as.matrix(pred)

  colnames(survmat) <- paste0("t=", times)

  as.data.frame(survmat)

  }



tune_selectcox <- function(formula, data, times,
                           rules = c("aic", "p"),
                           metrics = c("cindex", "ibs"),
                           folds = 5,
                           seed = 123,
                           refit_best = TRUE,
                           ...) {

  results <- purrr::map_dfr(rules, function(rule) {
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
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  results <- dplyr::arrange(results, dplyr::desc(!!sym(metrics[1])))

  if (refit_best) {
    best_rule <- results$rule[1]
    model <- fit_selectcox(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      rule = best_rule,
      ...
    )
    attr(model, "tuning_results") <- results
    return(model)
    }

  class(results) <- c("tuned_surv", class(results))
  return(results)
  }

res_selectcox <- tune_selectcox(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  rules = c("aic", "p"),
  metrics = c("cindex", "ibs", "ise"),
  folds = 3,
  refit_best = FALSE
  )

print(res_selectcox)

class(res_selectcox)

mod_selectcox <- tune_selectcox(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  rules = c("aic", "p"),
  metrics = c("cindex", "ibs", "ise"),
  folds = 3,
  refit_best = TRUE
  )

summary(mod_selectcox)


