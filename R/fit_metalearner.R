fit_metalearner_nnls_multi <- function(base_preds, time, status, times,
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

  structure(
    list(
      weights = weight_matrix,
      times = times,
      learners = learners,
      base_models = base_models,
      formula = formula,
      data = data
    ),
    class = "metalearner_nnls_multi",
    engine = "metalearner_nnls_multi"
  )
}

predict_metalearner_nnls_multi <- function(model, newdata, times) {
  stopifnot(inherits(model, "metalearner_nnls_multi"))

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

  colnames(out) <- requested_tnames
  as.data.frame(out)
}

plot_weights <- function(model) {
  stopifnot(inherits(model, "metalearner_nnls_multi"))
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

times <- c(1, 3, 5)
mod_aareg  <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)
mod_rpart  <- fit_rpart(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)
mod_ranger <- fit_ranger(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)
base_preds <- list(
  aareg  = predict_aareg(mod_aareg, sTRACE, times),
  rpart  = predict_rpart(mod_rpart, sTRACE, times),
  ranger = predict_ranger(mod_ranger, sTRACE, times)
)
meta_model <- fit_metalearner_nnls_multi(
  base_preds, sTRACE$time, sTRACE$status == 9, times,
  base_models = list(aareg = mod_aareg, rpart = mod_rpart, ranger = mod_ranger),
  formula = mod_aareg$formula, data = sTRACE
)

eval_all <- list(
  aareg  = evaluate_survlearner(mod_aareg,  metrics = c("cindex", "ibs", "iae", "ise"), times = times),
  rpart  = evaluate_survlearner(mod_rpart,  metrics = c("cindex", "ibs", "iae", "ise"), times = times),
  ranger = evaluate_survlearner(mod_ranger, metrics = c("cindex", "ibs", "iae", "ise"), times = times),
  meta   = evaluate_survlearner(meta_model, metrics = c("cindex", "ibs", "iae", "ise"), times = times)
)

do.call(rbind, eval_all)
plot_weights(meta_model)



#********



times <- c(1, 3, 5)
formula <- Surv(time, status == 9) ~ sex + diabetes + chf + vf

mod_aareg   <- fit_aareg(formula, data = sTRACE)
mod_rpart   <- fit_rpart(formula, data = sTRACE)
mod_ranger  <- fit_ranger(formula, data = sTRACE)
mod_cforest <- fit_cforest(formula, data = sTRACE)
mod_xgboost <- fit_xgboost(formula, data = sTRACE)

base_preds <- list(
  aareg   = predict_aareg(mod_aareg, sTRACE, times),
  rpart   = predict_rpart(mod_rpart, sTRACE, times),
  ranger  = predict_ranger(mod_ranger, sTRACE, times),
  cforest = predict_cforest(mod_cforest, sTRACE, times),
  xgboost = predict_xgboost(mod_xgboost, sTRACE, times))

base_models <- list(
  aareg = mod_aareg, rpart = mod_rpart, ranger = mod_ranger, cforest = mod_cforest,
  xgboost = mod_xgboost )

meta_model <- fit_metalearner_nnls_multi(
  base_preds, sTRACE$time, sTRACE$status == 9, times,
  base_models = base_models,
  formula = formula, data = sTRACE
)

eval_all <- lapply(base_models, function(mod) {
  evaluate_survlearner(mod, metrics = c("cindex", "ibs", "iae", "ise"), times = times)
})
eval_all$meta <- evaluate_survlearner(meta_model, metrics = c("cindex", "ibs", "iae", "ise"), times = times)

do.call(rbind, eval_all)
plot_weights(meta_model)

eval_all

library(dplyr)
tibble::tibble(
  learner = names(eval_all),
  cindex = lapply(eval_all, `[[`, "cindex")
) %>%
  tidyr::unnest(cols = cindex) %>%
  arrange(desc(cindex))
