fit_bart <- function(formula, data, K = 3, ...) {
  stopifnot(requireNamespace("BART", quietly = TRUE))

  mf <- model.frame(formula, data)
  y <- model.response(mf)

  if (!inherits(y, "Surv")) stop("Response must be a 'Surv' object.")
  if (attr(y, "type") != "right") stop("Only right-censored data supported.")

  time_vec   <- y[, "time"]
  status_vec <- y[, "status"]
  x <- model.matrix(formula, data = mf)[, -1, drop = FALSE]

  # silence verbose output
  invisible(
    capture.output({
      bart_fit <- BART::surv.bart(
        x.train = x,
        times = time_vec,
        delta = status_vec,
        K = K,
        printevery = Inf,  # optional
        ...
      )
    })
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

  # silently run prediction (suppresses C++ output)
  invisible(capture.output(
    pred_values <- predict(object$model, newdata = newx_expanded)$surv.test.mean
  ))

  prob_matrix <- matrix(pred_values, nrow = nrow(newx), ncol = K, byrow = TRUE)

  # align internal BART times to requested times via nearest neighbor match
  mapped_times <- sapply(times, function(t) which.min(abs(object$eval_times - t)))
  result <- prob_matrix[, mapped_times, drop = FALSE]
  colnames(result) <- paste0("t=", times)
  rownames(result) <- paste0("ID_", seq_len(nrow(newx)))

  return(as.data.frame(result))
}



library(survival)
library(BART)

# simulated example dataset
data(veteran)
veteran2 <- veteran
names(veteran2)[names(veteran2) == "time"] <- "Time"
names(veteran2)[names(veteran2) == "status"] <- "Event"

# fit and CV
mod_bart <- fit_bart(Surv(Time, Event) ~ age + karno + celltype, data = veteran2)

pred_bart <- predict_bart(mod_bart, newdata = veteran2[1:5, ], times = c(10, 30, 60))
print(round(pred_bart, 3))


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

# summary and plot
print(cv_results_bart)

# Summary and plot
print(cv_results_bart)
cv_summary(cv_results_bart)
cv_plot(cv_results_bart)

#---------------------- add tuner
tune_bart <- function(formula, data, times,
                      param_grid = expand.grid(
                        K = c(3, 5),
                        ntree = c(50, 100),
                        power = c(2, 2.5),
                        base = c(0.75, 0.95)
                      ),
                      metrics = c("cindex", "ibs"),
                      folds = 5, seed = 123,
                      refit_best = FALSE, ...) {

  results <- purrr::pmap_dfr(param_grid, function(K, ntree, power, base) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_bart,
      pred_fun = predict_bart,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      K = K, ntree = ntree, power = power, base = base,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(K = K, ntree = ntree, power = power, base = base) |>
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
    best_model <- fit_bart(
      formula = formula,
      data = data,
      K = best_row$K,
      ntree = best_row$ntree,
      power = best_row$power,
      base = best_row$base,
      ...
    )
    return(best_model)  # clean mlsurv_model
  }
}


param_grid = expand.grid(
  K = c(3),
  ntree = c(50),
  power = c(2),
  base = c(0.75, 0.95)
)

res <- tune_bart(
  formula = Surv(Time, Event) ~ age + karno + celltype,
  data = veteran2,
  param_grid = param_grid,
  times = c(10, 60),
  refit_best = FALSE
)

print(res)



mod <- tune_bart(
  formula = Surv(Time, Event) ~ age + karno + celltype,
  data = veteran2,
  param_grid = param_grid,
  times = c(10, 60),
  refit_best = TRUE
)

summary(mod)
predict_bart(mod, newdata = veteran2[1:5, ], times = c(30, 60))
