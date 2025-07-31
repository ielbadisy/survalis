fit_rfsrc <- function(formula, data,
                      ntree = 500,
                      mtry = NULL,
                      nodesize = 15,
                      ...) {
  stopifnot(requireNamespace("randomForestSRC", quietly = TRUE))

  model <- randomForestSRC::rfsrc(
    formula = formula,
    data = data,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    ...
  )

  structure(list(
    model = model,
    learner = "rfsrc",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "rfsrc")
}


predict_rfsrc <- function(object, newdata, times = NULL, ...) {
  stopifnot(object$learner == "rfsrc")
  stopifnot(requireNamespace("randomForestSRC", quietly = TRUE))

  pred <- randomForestSRC::predict.rfsrc(object$model, newdata = newdata)
  surv_mat <- pred$survival
  time_points <- pred$time.interest

  # interpolation to align predictions to `times`
  if (!is.null(times)) {
    idx <- sapply(times, function(t) which.min(abs(time_points - t)))
    surv_mat <- surv_mat[, idx, drop = FALSE]
    colnames(surv_mat) <- paste0("t=", times)
  } else {
    colnames(surv_mat) <- paste0("t=", time_points)
  }

  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(newdata)))
  return(as.data.frame(surv_mat))
}



library(survival)
library(randomForestSRC)

data(veteran, package = "survival")

mod_rsf <- fit_rfsrc(Surv(time, status) ~ age + celltype + karno, data = veteran, ntree = 200)

times <- c(100, 200, 300)
pred_probs <- predict_rfsrc(mod_rsf, newdata = veteran[1:5, ], times = times)

print(round(pred_probs, 3))

#------------- add tuner

tune_rfsrc <- function(formula, data, times,
                       param_grid = expand.grid(
                         ntree = c(200, 500),
                         mtry = c(1, 2, 3),
                         nodesize = c(5, 15)
                       ),
                       metrics = c("cindex", "ibs"),
                       folds = 5,
                       seed = 123,
                       refit_best = FALSE,
                       ...) {

  results <- purrr::pmap_dfr(param_grid, function(ntree, mtry, nodesize) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_rfsrc,
      pred_fun = predict_rfsrc,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))  # like other learners
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_rfsrc(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      ntree = best_row$ntree,
      mtry = best_row$mtry,
      nodesize = best_row$nodesize,
      ...
    )
    attr(best_model, "tuning_results") <- results  # optional metadata
    return(best_model)
  }
}


mod_rfsrc_best <- tune_rfsrc(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    ntree = c(200, 500),
    mtry = c(1, 2, 3),
    nodesize = c(5, 15)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

summary(mod_rfsrc_best)
