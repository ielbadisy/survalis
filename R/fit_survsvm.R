fit_survsvm <- function(formula, data,
                            type = "regression",
                            gamma.mu = 0.1,
                            opt.meth = "quadprog",
                            kernel = "add_kernel",
                            diff.meth = NULL,
                            ...) {

  stopifnot(requireNamespace("survivalsvm", quietly = TRUE))

  model <- survivalsvm::survivalsvm(
    formula = formula,
    data = data,
    type = type,
    gamma.mu = gamma.mu,
    opt.meth = opt.meth,
    kernel = kernel,
    diff.meth = diff.meth,
    ...
  )

  structure(list(
    model = model,
    learner = "survsvm",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "survivalsvm")
}





predict_survsvm <- function(object, newdata, times, dist = "exp", shape = 1) {

  if (!is.null(object$learner) && object$learner != "survsvm") {
    warning("Object passed to predict_survsvm() may not come from fit_survsvm().")
  }

  stopifnot(requireNamespace("survivalsvm", quietly = TRUE))

  newdata <- as.data.frame(newdata)
  pred <- predict(object$model, newdata = newdata)
  predicted_times <- as.numeric(pred$predicted)

  dist <- match.arg(dist, choices = c("exp", "weibull"))

  survmat <- switch(dist,
                  "exp" = sapply(times, function(t) exp(-t / predicted_times)),
                  "weibull" = sapply(times, function(t) exp(- (t / predicted_times)^shape))
  )

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}

library(survival)
library(survivalsvm)
library(ggplot2)

mod_svm <- fit_survsvm(Surv(time, status) ~ age + celltype + karno,
                           data = veteran,
                           type = "regression",
                           gamma.mu = 0.1,
                           kernel = "lin_kernel")

times <- c(100, 300, 500)
predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "exp")
predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = times, dist = "weibull", shape = 1.5)

#-------------------

cv_results_svm <- cv_survlearner(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  fit_fun = fit_survsvm,
  pred_fun = predict_survsvm,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42,
  gamma.mu = 0.1,
  kernel = "lin_kernel")

print(cv_results_svm)
cv_summary(cv_results_svm)
cv_plot(cv_results_svm)


#------- add tuner

tune_survsvm <- function(formula, data, times,
                             metrics = "cindex",
                             param_grid,
                             folds = 5, seed = 42,
                             refit_best = FALSE,
                             dist = "exp", shape = 1) {
  stopifnot(is.list(param_grid))
  param_df <- tidyr::crossing(!!!param_grid)

  # Fixed safe defaults
  fixed_type <- "regression"
  fixed_opt_meth <- "quadprog"
  fixed_diff_meth <- NULL

  results <- purrr::pmap_dfr(param_df, function(gamma.mu, kernel) {
    set.seed(seed)
    tryCatch({
      res_cv <- cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_survsvm,
        pred_fun = function(...) predict_survsvm(..., dist = dist, shape = shape),
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        gamma.mu = gamma.mu,
        kernel = kernel,
        type = fixed_type,
        opt.meth = fixed_opt_meth,
        diff.meth = fixed_diff_meth
      )

      summary <- cv_summary(res_cv)

      tibble::tibble(
        gamma.mu = gamma.mu,
        kernel = kernel
      ) |>
        dplyr::bind_cols(
          tidyr::pivot_wider(summary[, c("metric", "mean")],
                             names_from = metric, values_from = mean)
        )
    }, error = function(e) {
      invisible(cat(glue::glue("Skipping (gamma.mu={gamma.mu}, kernel={kernel}): {e$message}\n"),
                    file = "errors.log", append = TRUE)) ## keep a trace without showing it on the console
      NULL
    })
  })

  if (nrow(results) == 0) {
    warning("All tuning combinations failed.")
    return(NULL)
  }

  results <- dplyr::arrange(results, dplyr::desc(.data[[metrics[1]]]))

  if (refit_best) {
    best <- results[1, ]
    model <- fit_survsvm(
      formula = formula,
      data = data,
      gamma.mu = best$gamma.mu,
      kernel = best$kernel,
      type = fixed_type,
      opt.meth = fixed_opt_meth,
      diff.meth = fixed_diff_meth
    )
    attr(model, "tuning_results") <- results
    class(model) <- c("tune_surv", class(model))
    return(model)
  }

  class(results) <- c("tune_surv", class(results))
  return(results)
}




grid <- list(
  gamma.mu = c(0.01, 0.1),
  kernel = c("lin_kernel", "add_kernel")
)


res_svm <- tune_survsvm(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs"),
  param_grid = grid,
  folds = 3,
  refit_best = TRUE
)

summary(res_svm)



res_svm <- tune_survsvm(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs"),
  param_grid = grid,
  folds = 3,
  refit_best = FALSE
)

res_svm
