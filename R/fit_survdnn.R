fit_survdnn <- function(formula, data,
                        loss = "cox",
                        hidden = c(32L, 32L, 16L),
                        activation = "relu",
                        lr = 1e-4,
                        epochs = 300L,
                        verbose = FALSE,
                        ...) {

  stopifnot(requireNamespace("survdnn", quietly = TRUE))
  stopifnot(requireNamespace("torch", quietly = TRUE))

  model <- survdnn::survdnn(
    formula = formula,
    data = data,
    loss = loss,
    hidden = hidden,
    activation = activation,
    lr = lr,
    epochs = epochs,
    verbose = verbose
  )

  time_status <- all.vars(formula[[2]])

  structure(list(
    model = model,
    learner = "survdnn",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2]
  ), class = "mlsurv_model", engine = "survdnn")
  }





predict_survdnn <- function(object, newdata, times, ...) {


  if (!is.null(object$learner) && object$learner != "survdnn") {
    warning("Object passed to predict_survdnn() may not come from fit_survdnn().")
  }

  predict(object$model, newdata = newdata, type = "survival", times = times)
}


library(survdnn)
library(survival)
data(veteran)

mod <- fit_survdnn(Surv(time, status) ~ age + karno + celltype,
                   data = veteran, loss = "cox", epochs = 50, verbose = FALSE)

pred <- predict_survdnn(mod, newdata = veteran[1:5, ], times = c(30, 90, 180))
print(pred)



tune_survdnn <- function(formula, data, times,
                         param_grid = list(
                           hidden = list(c(32, 16), c(64, 32, 16)),
                           lr = c(1e-3, 5e-4),
                           activation = c("relu", "gelu"),
                           epochs = c(100, 200),
                           loss = c("cox", "aft")
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 3,
                         seed = 42,
                         refit_best = FALSE,
                         ...) {

  stopifnot(!missing(formula), !missing(data), !missing(times), !missing(param_grid))
  stopifnot(is.list(param_grid))

  param_df <- tidyr::crossing(!!!param_grid)

  results <- purrr::pmap_dfr(param_df, function(hidden, lr, activation, epochs, loss) {
    set.seed(seed)
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_survdnn,
      pred_fun = predict_survdnn,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      hidden = hidden,
      lr = lr,
      activation = activation,
      epochs = epochs,
      loss = loss,
      verbose = FALSE,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      hidden = list(hidden),
      lr = lr,
      activation = activation,
      epochs = epochs,
      loss = loss
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric,
                           values_from = mean)
      )
  })

  results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_survdnn(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      hidden = best_row$hidden[[1]],
      lr = best_row$lr,
      activation = best_row$activation,
      epochs = best_row$epochs,
      loss = best_row$loss,
      verbose = FALSE,
      ...
    )
    attr(best_model, "tuning_results") <- results
    return(best_model)
  }
}



library(survival)
data(veteran)

grid <- list(
  hidden = list(c(16), c(32, 16)),
  lr = c(1e-4, 5e-4),
  activation = c("relu", "tanh"),
  epochs = c(300),
  loss = c("cox", "coxtime")

  )



mod <- tune_survdnn(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(90),
  metrics = c("cindex", "ibs"),
  param_grid = grid,
  refit_best = TRUE
)




summary(mod)