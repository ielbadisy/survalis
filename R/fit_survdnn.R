
fit_survdnn <- function(formula, data,
                        loss = "cox",
                        hidden = c(32L, 32L, 16L),
                        activation = "relu",
                        lr = 1e-4,
                        epochs = 300L,
                        verbose = TRUE,
                        ...) {

  stopifnot(requireNamespace("survdnn", quietly = TRUE))
  stopifnot(requireNamespace("torch", quietly = TRUE))


  model <- survdnn(
    formula = formula,
    data = data,
    loss = loss,
    hidden = hidden,
    activation = activation,
    lr = lr,
    epochs = epochs,
    verbose = verbose
  )

  result <- list(
    model = model,
    formula = formula,
    data = data
  )

  attr(result, "engine") <- "survdnn"
  class(result) <- "mlsurv_model"
  return(result)
}



#' Predict Survival Probabilities from survdnn Model
#'
#' @param model A model fitted with `fit_survdnn`
#' @param newdata A data frame with covariates
#' @param times A numeric vector of prediction time points
#' @param ... Not used
#'
#' @return A data frame (nrow = nrow(newdata), ncol = length(times)) of survival probs
#' @export
predict_survdnn <- function(model, newdata, times, ...) {
  predict(model$model, newdata = newdata, type = "survival", times = times)
}


library(survdnn)
library(survival)
data(veteran)

mod <- fit_survdnn(Surv(time, status) ~ age + karno + celltype,
                   data = veteran, loss = "cox", epochs = 50, verbose = FALSE)

pred <- predict_survdnn(mod, newdata = veteran[1:5, ], times = c(30, 90, 180))
print(pred)

# Evaluate C-index, IBS, IAE, ISE (multi-time OK)
evaluate_survlearner(
  model = mod,
  metrics = c("cindex", "ibs", "iae", "ise"),
  times = c(30, 90, 180)
)

# Evaluate Brier (must use only one time point)
evaluate_survlearner(
  model = mod,
  metrics = "brier",
  times = 90
)




tune_survdnn <- function(formula, data, times, metrics = "cindex",
                         param_grid, cv = 3, seed = 42) {
  stopifnot(!missing(formula), !missing(data), !missing(times), !missing(param_grid))
  stopifnot(is.list(param_grid), "hidden" %in% names(param_grid), "lr" %in% names(param_grid))

  param_df <- tidyr::crossing(!!!param_grid)

  results <- purrr::pmap_dfr(param_df, function(hidden, lr, activation, epochs, loss = "cox") {
    set.seed(seed)
    mod <- fit_survdnn(
      formula = formula,
      data = data,
      hidden = hidden,
      lr = lr,
      activation = activation,
      epochs = epochs,
      loss = loss,
      verbose = FALSE
    )

    metrics_tbl <- evaluate_survlearner(mod, metrics = metrics, times = times)

    tibble::tibble(
      hidden = list(hidden),
      lr = lr,
      activation = activation,
      epochs = epochs,
      loss = loss,
      metric = metrics_tbl$metric,
      value = metrics_tbl$value
    )
  })

  results %>% dplyr::arrange(metric, dplyr::desc(value))
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

tune_survdnn(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(90),
  metrics = c("cindex", "ibs"),
  param_grid = grid
)




