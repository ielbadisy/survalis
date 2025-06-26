library(data.table)
fit_bart <- function(formula, data, K = 3, ...) {
  stopifnot(requireNamespace("BART", quietly = TRUE))

  mf <- model.frame(formula, data)
  y <- model.response(mf)

  if (!inherits(y, "Surv")) stop("Response must be a 'Surv' object.")
  if (attr(y, "type") != "right") stop("Only right-censored data supported.")

  time_vec   <- y[, "time"]
  status_vec <- y[, "status"]

  x <- model.matrix(formula, data = mf)[, -1, drop = FALSE]  # drop intercept

  ## need to update re-implment the new version of this algorithm
  bart_fit <- BART::surv.bart(
    x.train = x,
    times = time_vec,
    delta = status_vec,
    K = K,
    printevery = 10000, # or Inf
    ...
  )

  structure(
    list(
      model = bart_fit,
      learner = "bart",
      formula = formula,
      data = data,
      eval_times = bart_fit$times
    ),
    class = "mlsurv_model",
    engine = "bart"
  )
}


predict_bart <- function(object, newdata, times, ...) {
  if (object$learner != "bart") {
    stop("This model is not a BART model.")
  }

  if (missing(times)) stop("Please provide `times` for prediction.")
  stopifnot(requireNamespace("BART", quietly = TRUE))

  newx <- model.matrix(object$formula, data = newdata)[, -1, drop = FALSE]
  K <- length(object$eval_times)

  newx_expanded <- cbind(
    t = object$eval_times,
    newx[rep(seq_len(nrow(newx)), each = K), ]
  )

  pred <- predict(object$model, newdata = newx_expanded)$surv.test.mean
  prob_matrix <- matrix(pred, nrow = nrow(newx), ncol = K, byrow = TRUE)

  # align internal BART times to requested times via nearest neighbor match
  mapped_times <- sapply(times, function(t) which.min(abs(object$eval_times - t)))
  closest_times <- object$eval_times[mapped_times]
  result <- prob_matrix[, mapped_times, drop = FALSE]
  colnames(result) <- paste0("t=", times)  # preserve user intent

  rownames(result) <- paste0("ID_", seq_len(nrow(newx)))

  return(as.data.frame(result))
}



library(survival)
library(BART)

# Simulated example dataset
data(veteran)
veteran2 <- veteran
names(veteran2)[names(veteran2) == "time"] <- "Time"
names(veteran2)[names(veteran2) == "status"] <- "Event"

# Fit and CV
mod_bart <- fit_bart(Surv(Time, Event) ~ age + karno + celltype, data = veteran2)





cv_results_bart <- cv_survlearner(
  formula = Surv(Time, Event) ~ age + karno + celltype,
  data = veteran2,
  fit_fun = fit_bart,
  pred_fun = predict_bart,
  times = c(5, 10, 40),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42
)

# Summary and plot
print(cv_results_bart)
cv_summary(cv_results_bart)
cv_plot(cv_results_bart)

#---------------------- add tuner


#' Tune Bayesian Additive Regression Trees for Survival (BART)
#'
#' @param formula A survival formula
#' @param data A data.frame
#' @param times Time points to evaluate survival probabilities
#' @param param_grid A data.frame of tuning grid (K, ntree, etc.)
#' @param metrics Vector of evaluation metrics (e.g., "cindex", "ibs")
#' @param folds Number of CV folds
#' @param seed Seed for reproducibility
#' @param ... Extra arguments passed to `fit_bart()`
#'
#' @return A tibble with average performance across CV folds
#' @export
tune_bart <- function(formula, data, times,
                      param_grid = expand.grid(
                        K = c(3, 5),
                        ntree = c(50, 100),
                        power = c(2, 2.5),
                        base = c(0.75, 0.95)
                      ),
                      metrics = c("cindex", "ibs"),
                      folds = 5, seed = 123, ...) {

  purrr::pmap_dfr(param_grid, function(K, ntree, power, base) {

    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_bart,
      pred_fun = predict_bart,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      K = K,
      ntree = ntree,
      power = power,
      base = base,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(K = K, ntree = ntree, power = power, base = base) |>
      bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric,
                           values_from = mean)
      )
  }) |> arrange(dplyr::across(any_of(metrics[1])))
}



res_bart <- tune_bart(
  formula = Surv(Time, Event) ~ age + karno + celltype,
  data = veteran2,
  times = c(10, 60),
  param_grid = expand.grid(
    K = c(5, 10),
    ntree = c(10),
    power = c(3),
    base = c(0.95)
  ),
  folds = 3,
  metrics = c("cindex", "ibs")
)

print(res_bart)

