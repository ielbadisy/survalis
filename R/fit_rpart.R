fit_rpart <- function(formula, data,
                      minsplit = 20,
                      minbucket = round(minsplit / 3),
                      cp = 0.01,
                      maxcompete = 4,
                      maxsurrogate = 5,
                      usesurrogate = 2,
                      xval = 10,
                      surrogatestyle = 0,
                      maxdepth = 30) {
  stopifnot(requireNamespace("rpart", quietly = TRUE))

  method <- "exp"  # fixed to survival tree with exponential splitting

  model <- rpart::rpart(
    formula,
    data = data,
    method = method,
    control = list(
      minsplit = minsplit, minbucket = minbucket, cp = cp,
      maxcompete = maxcompete, maxsurrogate = maxsurrogate,
      usesurrogate = usesurrogate, xval = xval,
      surrogatestyle = surrogatestyle, maxdepth = maxdepth
    ),
    model = TRUE, x = TRUE, y = TRUE
  )

  structure(list(
    model = model,
    learner = "rpart",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2]
  ), class = "mlsurv_model", engine = "rpart")
}

predict_rpart <- function(object, newdata, times, ...) {

  if (!is.null(object$learner) && object$learner != "rpart") {
    warning("Object passed to predict_rpart() may not come from fit_rpart().")
  }

  stopifnot(object$learner == "rpart")
  stopifnot(requireNamespace("partykit", quietly = TRUE))

  model_party <- partykit::as.party(object$model)

  # predict expected survival times from the tree
  predicted_times <- stats::predict(model_party, newdata = newdata, type = "response")

  # convert predicted survival times to survival probabilities using exponential assumption
  survmat <- sapply(times, function(t) {
    exp(-t / predicted_times)
  })

  if (is.vector(survmat)) survmat <- matrix(survmat, ncol = length(times))

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}



library(survival)
library(rpart)
library(partykit)

mod_rpart <- fit_rpart(Surv(time, status) ~ age + karno + celltype, data = veteran)

pred_rpart <- predict_rpart(mod_rpart, newdata = veteran[1:5, ], times = c(100, 200, 300, 400, 500))


## predict_rpart has a weird behavior: predicted probs are repetitives -> maybe this came from the fact that this learner predict one single prob for all individuals in the same node

pred_rpart



cv_results_rpart <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  fit_fun = fit_rpart,
  pred_fun = predict_rpart,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42
)

print(cv_results_rpart)
cv_summary(cv_results_rpart)
cv_plot(cv_results_rpart)

#------------- add tuner

tune_rpart <- function(formula, data, times,
                       param_grid = expand.grid(
                         minsplit = c(10, 20),
                         cp = c(0.001, 0.01),
                         maxdepth = c(10, 30)
                       ),
                       metrics = c("cindex", "ibs"),
                       folds = 5,
                       seed = 123,
                       refit_best = TRUE,
                       ...) {
  res <- purrr::pmap_dfr(param_grid, function(minsplit, cp, maxdepth) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_rpart,
      pred_fun = predict_rpart,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      minsplit = minsplit,
      cp = cp,
      maxdepth = maxdepth,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(
      minsplit = minsplit,
      cp = cp,
      maxdepth = maxdepth
    ) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric, values_from = mean)
      )
  })

  res <- dplyr::arrange(res, dplyr::desc(!!sym(metrics[1])))

  if (refit_best) {
    best_params <- res[1, ]
    model <- fit_rpart(
      formula = formula,
      data = tidyr::drop_na(data, all.vars(formula)),
      minsplit = best_params$minsplit,
      cp = best_params$cp,
      maxdepth = best_params$maxdepth,
      ...
    )
    attr(model, "tuning_results") <- res
    class(model) <- c("tune_surv", class(model))
    return(model)
  }

  return(res)
}



# Load required packages
library(survival)
library(rpart)
library(partykit)
library(dplyr)

# Source or define fit_rpart() and predict_rpart() beforehand
# Assuming your definitions are already loaded

# Load the dataset


veteran <- survival::veteran
# Set up a formula and time points
form <- Surv(time, status) ~ age + karno + celltype
times <-1:100

# 1. Test tuning WITHOUT refitting
res_rpart <- tune_rpart(
  formula = form,
  data = veteran,
  times = times,
  param_grid = expand.grid(
    minsplit = c(10, 20),
    cp = c(0.001, 0.01),
    maxdepth = c(10, 30)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  seed = 42,
  refit_best = FALSE
)

# Should print a tibble of cross-validated scores
print(res_rpart)
class(res_rpart)

# 2. Test tuning WITH refitting of the best model
mod_rpart_best <- tune_rpart(
  formula = form,
  data = veteran,
  times = times,
  param_grid = expand.grid(
    minsplit = c(10, 20),
    cp = c(0.001, 0.01),
    maxdepth = c(10, 30)
  ),
  metrics = c("cindex", "ibs"),
  folds = 3,
  seed = 42,
  refit_best = TRUE
)

# Should return a fitted model object
print(mod_rpart_best)
summary(mod_rpart_best)

# 3. Predict survival probabilities using the best model
pred_rpart <- predict_rpart(mod_rpart_best, newdata = veteran[1:5, ], times = times)
print(round(pred_rpart, 3))

