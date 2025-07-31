fit_cforest <- function(formula, data,
                        teststat = "quad",
                        testtype = "Univariate",
                        mincriterion = 0,
                        ntree = 500,
                        mtry = 5,
                        replace = TRUE,
                        fraction = 0.632,
                        ...) {
  stopifnot(requireNamespace("party", quietly = TRUE))

  model <- party::cforest(
    formula = formula,
    data = data,
    controls = party::cforest_control(
      teststat = teststat,
      testtype = testtype,
      mincriterion = mincriterion,
      ntree = ntree,
      mtry = mtry,
      replace = replace,
      fraction = fraction,
      ...
    )
  )

  structure(
    list(
      model = model,
      learner = "cforest",
      formula = formula,
      data = data
    ),
    class = "mlsurv_model",
    engine = "cforest"
  )
}


predict_cforest <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "cforest")
  stopifnot(requireNamespace("party", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  # get survival predictions from cforest
  survlist <- suppressMessages(suppressWarnings(
    predict(object$model, newdata = newdata, type = "prob")
  ))

  # interpolate at requested times
  prob_matrix <- t(sapply(survlist, function(sfit) {
    s_time <- sfit$time
    s_surv <- sfit$surv

    if (is.null(s_time) || is.null(s_surv) || length(s_surv) == 0) {
      return(rep(NA, length(times)))
    }

    # Linear interpolation with rule = 2 for extrapolation
    approx(x = s_time, y = s_surv, xout = times, method = "linear", rule = 2)$y
  }))

  # label columns and rows as in other learners
  colnames(prob_matrix) <- paste0("t=", times)
  rownames(prob_matrix) <- paste0("ID_", seq_len(nrow(newdata)))

  return(as.data.frame(prob_matrix))
}


library(survival)
data(veteran)

mod_cforest <- fit_cforest(Surv(time, status) ~ age + celltype + karno, data = veteran)
predict_cforest(mod_cforest, newdata = veteran[1:5, ], times = c(100, 200, 300))



# --- Evaluate with CV ---
cv_results_cforest <- cv_survlearner(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  fit_fun = fit_cforest,
  pred_fun = predict_cforest,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 10,
  seed = 42
)

# --- Print Results ---
print(cv_results_cforest)



#---------------- tune_cforest()

tune_cforest <- function(formula, data, times,
                         param_grid = expand.grid(
                           ntree = c(100, 300),
                           mtry = c(2, 3),
                           mincriterion = c(0, 0.95),
                           fraction = c(0.5, 0.632)
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5, seed = 123,
                         refit_best = FALSE, ...) {

  results <- purrr::pmap_dfr(param_grid, function(ntree, mtry, mincriterion, fraction) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_cforest,
      pred_fun = predict_cforest,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ntree = ntree,
      mtry = mtry,
      mincriterion = mincriterion,
      fraction = fraction,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(ntree = ntree,
                   mtry = mtry,
                   mincriterion = mincriterion,
                   fraction = fraction) |>
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
    best_model <- fit_cforest(
      formula = formula,
      data = data,
      ntree = best_row$ntree,
      mtry = best_row$mtry,
      mincriterion = best_row$mincriterion,
      fraction = best_row$fraction,
      ...
    )
    return(best_model)
  }
}


# Define parameter grid
param_grid <- expand.grid(
  ntree = c(100),
  mtry = c(2),
  mincriterion = c(0),
  fraction = c(0.632)
)

# Return CV results only
res_cforest <- tune_cforest(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = param_grid,
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = FALSE
)

print(res_cforest)

# Return best fitted model directly (mlsurv_model)
mod_cforest <- tune_cforest(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  times = c(100, 200, 300),
  param_grid = param_grid,
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

summary(mod_cforest)
predict_cforest(mod_cforest, newdata = veteran[1:5, ], times = c(100, 200, 300))



