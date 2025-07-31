
fit_ranger <- function(formula, data, ...) {
  stopifnot(requireNamespace("ranger", quietly = TRUE))

  model <- ranger::ranger(
    formula = formula,
    data = data,
    ...,
    respect.unordered.factors = "order",
    case.weights = NULL
  )

  structure(
    list(
      model = model,
      learner = "ranger",
      formula = formula,
      data = data
    ),
    class = "mlsurv_model",
    engine = "ranger"
  )
}


predict_ranger <- function(model, newdata, times) {
  requireNamespace("ranger")

  pred_obj <- predict(model$model, data = newdata)
  pred <- pred_obj$survival
  model_times <- model$model$unique.death.times

  if (is.null(pred)) stop("Prediction returned NULL survival matrix.")
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)

  surv_mat <- matrix(NA_real_, nrow = nrow(pred), ncol = length(times))
  for (i in seq_len(nrow(pred))) {
    y <- pred[i, ]
    surv_mat[i, ] <- stats::approx(
      x = model_times, y = y, xout = times,
      method = "linear", rule = 2, ties = "ordered"
    )$y
  }

  colnames(surv_mat) <- paste0("t=", times)
  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(surv_mat)))
  as.data.frame(surv_mat)
}



mod <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
summary(mod)

tune_ranger <- function(formula, data, times,
                        param_grid = expand.grid(
                          num.trees = c(100, 300),
                          mtry = c(1, 2, 3),
                          min.node.size = c(3, 5)
                        ),
                        metrics = c("cindex", "ibs"),
                        folds = 5, seed = 123,
                        refit_best = FALSE,
                        ...) {

  # cv loop
  res <- purrr::pmap_dfr(param_grid, function(num.trees, mtry, min.node.size) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_ranger,
      pred_fun = predict_ranger,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      num.trees = num.trees,
      mtry = mtry,
      min.node.size = min.node.size,
      ...
    )

    summary <- cv_summary(cv_results)
    tibble::tibble(num.trees, mtry, min.node.size) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  res <- dplyr::arrange(res, dplyr::across(dplyr::any_of(metrics[1])))

  if (refit_best) {
    best <- res[1, ]
    mod <- fit_ranger(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      num.trees = best$num.trees,
      mtry = best$mtry,
      min.node.size = best$min.node.size,
      ...
    )
    attr(mod, "tuning_results") <- res
    class(mod) <- c("tune_surv", class(mod))  # this line enables summary()
    return(mod)
  }

  return(res)
}



mod_ranger_best <- tune_ranger(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    num.trees = c(100, 300),
    mtry = c(1, 2),
    min.node.size = c(3, 5)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

# This now works and gives structured summary output:
summary(mod_ranger_best)




mod_ranger <- tune_ranger(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = expand.grid(
    num.trees = c(100, 300),
    mtry = c(1, 2),
    min.node.size = c(3, 5)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = FALSE
)

mod_ranger


