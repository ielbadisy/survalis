fit_metalearner <- function(base_preds, time, status, times,
                            base_models, formula, data) {
  stopifnot(is.list(base_preds), is.list(base_models), is.data.frame(data))

  base_preds <- lapply(base_preds, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })

  learners <- names(base_preds)
  T <- length(times)
  weight_matrix <- matrix(NA_real_, nrow = length(learners), ncol = T,
                          dimnames = list(learners, paste0("t=", times)))

  for (j in seq_along(times)) {
    t_star <- times[j]
    Y <- as.integer(time > t_star)
    X <- sapply(base_preds, function(S) S[, j])
    fit <- nnls::nnls(as.matrix(X), Y)
    w <- coef(fit); w <- w / sum(w)
    weight_matrix[, j] <- w
  }

  structure(list(
    weights = weight_matrix,
    base_models = base_models,
    base_preds = base_preds,
    formula = formula,
    data = data,
    time = time,
    status = status,
    learners = learners
  ), class = "metalearner", learner = "metalearner", engine = "survalis")
}


predict_metalearner <- function(model, newdata, times) {
  stopifnot(inherits(model, "metalearner"))

  learners <- model$learners
  base_models <- model$base_models
  requested_tnames <- paste0("t=", times)
  W <- model$weights[, requested_tnames, drop = FALSE]

  base_preds_test <- lapply(learners, function(learner) {
    pred_fun <- get(paste0("predict_", learner), mode = "function")
    pred <- pred_fun(base_models[[learner]], newdata = newdata, times = times)
    pred <- as.matrix(pred)
    colnames(pred) <- requested_tnames
    pred
  })
  names(base_preds_test) <- learners

  n <- nrow(base_preds_test[[1]])
  T <- length(times)
  pred_array <- array(NA_real_, dim = c(n, length(learners), T))
  for (k in seq_along(learners)) {
    pred_array[, k, ] <- base_preds_test[[k]]
  }

  out <- matrix(NA_real_, nrow = n, ncol = T)
  for (j in seq_len(T)) {
    out[, j] <- pred_array[, , j] %*% W[, j]
  }

  #colnames(out) <- requested_tnames
  as.data.frame(out)
}

plot_metalearner_weights <- function(model) {
  stopifnot(inherits(model, "metalearner"))
  library(ggplot2); library(tidyr); library(dplyr)

  W <- as.data.frame(model$weights)
  W$learner <- rownames(model$weights)

  W_long <- W |>
    pivot_longer(-learner, names_to = "time", values_to = "weight") |>
    mutate(time = as.numeric(sub("t=", "", time)))

  ggplot(W_long, aes(x = time, y = weight, color = learner)) +
    geom_line(size = 1.2) +
    labs(x = "Time", y = "Weight", title = "NNLS Stacking Weights Over Time") +
    theme_minimal(base_size = 14)
}



## TEST -------

library(survival)
data(veteran)

formula <- Surv(time, status) ~ age + karno + trt + diagtime + celltype
times <- c(30, 60, 90, 120, 180)

mod_aareg   <- fit_aareg(formula, data = veteran)
mod_rpart   <- fit_rpart(formula, data = veteran)
mod_ranger  <- fit_ranger(formula, data = veteran)
mod_cforest <- fit_cforest(formula, data = veteran)
mod_xgboost <- fit_xgboost(formula, data = veteran)

base_preds <- list(
  aareg   = predict_aareg(mod_aareg, veteran, times),
  rpart   = predict_rpart(mod_rpart, veteran, times),
  ranger  = predict_ranger(mod_ranger, veteran, times),
  cforest = predict_cforest(mod_cforest, veteran, times),
  xgboost = predict_xgboost(mod_xgboost, veteran, times)
)

base_models <- list(
  aareg   = mod_aareg,
  rpart   = mod_rpart,
  ranger  = mod_ranger,
  cforest = mod_cforest,
  xgboost = mod_xgboost
)

meta_model <- fit_metalearner(
  base_preds = base_preds,
  time = veteran$time,
  status = veteran$status,
  times = times,
  base_models = base_models,
  formula = formula,
  data = veteran
)

# Prediction
pred_meta <- predict_metalearner(meta_model, newdata = veteran[1:5, ], times = times)
print(round(pred_meta, 3))

# Plot
plot_metalearner_weights(meta_model)



eval_meta <- evaluate_survlearner(meta_model, metrics = c("cindex", "ibs", "iae", "ise"), times = times)


eval_meta
