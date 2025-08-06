fit_survmetalearner <- function(base_preds, time, status, times,
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
    learner = "survmetalearner",         # for predict_* lookup
    learners = learners
  ), class = c("mlsurv_model", "survmetalearner"), engine = "survalis")  # dual class

}


predict_survmetalearner <- function(model, newdata, times) {
  stopifnot(inherits(model, "survmetalearner"))

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

  survmat <- matrix(NA_real_, nrow = n, ncol = T)
  for (j in seq_len(T)) {
    survmat[, j] <- pred_array[, , j] %*% W[, j]
  }

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}

plot_survmetalearner_weights <- function(model) {
  stopifnot(inherits(model, "survmetalearner"))
  library(ggplot2); library(tidyr); library(dplyr)

  W <- as.data.frame(model$weights)
  W$learner <- rownames(model$weights)

  W_long <- W |>
    pivot_longer(-learner, names_to = "time", values_to = "weight") |>
    mutate(time = as.numeric(sub("t=", "", time)))

  ggplot(W_long, aes(x = time, y = weight, color = learner)) +
    geom_line(size = 1.2) +
    labs(x = "Time", y = "Weight", title = "NNLS stacking weights over time") +
    theme_minimal(base_size = 14)
}



## TEST -------

library(survival)
data(veteran)

formula <- Surv(time, status) ~ age + karno + trt + diagtime + celltype
times <- c(200, 500, 800)

mod_rpart   <- fit_rpart(formula, data = veteran)
mod_ranger  <- fit_ranger(formula, data = veteran)
mod_cforest <- fit_cforest(formula, data = veteran)

base_preds <- list(
  rpart   = predict_rpart(mod_rpart, veteran, times),
  ranger  = predict_ranger(mod_ranger, veteran, times),
  cforest = predict_cforest(mod_cforest, veteran, times)
)

base_models <- list(
  rpart   = mod_rpart,
  ranger  = mod_ranger,
  cforest = mod_cforest
)

meta_model <- fit_survmetalearner(
  base_preds = base_preds,
  time = veteran$time,
  status = veteran$status,
  times = times,
  base_models = base_models,
  formula = formula,
  data = veteran
)

# Prediction
pred_meta <- predict_survmetalearner(meta_model, newdata = veteran[1:5, ], times = times)
print(round(pred_meta, 3))

# Plot
plot_survmetalearner_weights(meta_model)



summary(meta_model)



## Score


score_survmodel(meta_model, times = c(200, 500, 800), metrics = c("cindex", "ibs"))



## CV


cv_survmetalearner <- function(formula, data, times,
                           base_models, base_preds,
                           folds = 5,
                           metrics = c("cindex", "ibs"),
                           seed = 123, verbose = TRUE) {
  stopifnot(is.list(base_models), is.list(base_preds))
  stopifnot(inherits(formula, "formula"), is.data.frame(data))

  # Convert all base_preds to matrices (if needed)
  base_preds <- lapply(base_preds, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })

  # Validate shape
  n_obs <- nrow(data)
  for (learner in names(base_preds)) {
    if (nrow(base_preds[[learner]]) != n_obs) {
      stop(sprintf("Number of rows in base_preds[['%s']] does not match data.", learner))
    }
  }

  set.seed(seed)
  rsplits <- rsample::vfold_cv(data, v = folds)

  results <- purrr::map_dfr(seq_along(rsplits$splits), function(i) {
    split <- rsplits$splits[[i]]
    train <- rsample::analysis(split)
    test  <- rsample::assessment(split)

    if (verbose) cat("Processing fold", i, "...\n")

    # subset base_preds for training rows
    idx_train <- as.integer(rownames(data)) %in% as.integer(rownames(train))
    base_preds_train <- lapply(base_preds, function(mat) mat[idx_train, , drop = FALSE])

    # extract time/status from formula
    tf <- terms(formula, data = train)
    outcome <- attr(tf, "variables")[[2]]
    time_col <- as.character(outcome[[2]])
    status_expr <- outcome[[3]]

    if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
      status_col <- as.character(status_expr[[2]])
      event_value <- eval(status_expr[[3]], train)
      status_vector <- as.integer(train[[status_col]] == event_value)
    } else {
      status_col <- as.character(status_expr)
      event_value <- 1
      status_vector <- train[[status_col]]
    }

    # fit survmetalearner on training fold
    meta_fit <- fit_survmetalearner(
      base_preds = base_preds_train,
      time = train[[time_col]],
      status = status_vector,
      times = times,
      base_models = base_models,
      formula = formula,
      data = train
    )

    # predict on test
    preds_test <- predict_survmetalearner(meta_fit, newdata = test, times = times)

    # extract test survival object
    time_test <- test[[time_col]]
    status_test <- if (is.call(status_expr)) {
      as.integer(test[[as.character(status_expr[[2]])]] == eval(status_expr[[3]], test))
    } else {
      test[[as.character(status_expr)]]
    }
    surv_test <- survival::Surv(time_test, status_test)

    # score metrics
    tibble::tibble(
      fold = i,
      metric = metrics
    ) |>
      dplyr::mutate(value = purrr::map(metric, function(metric) {
        switch(metric,
               "cindex" = cindex_survmat(surv_test, predicted = preds_test, t_star = max(times)),
               "brier"  = if (length(times) != 1) stop("Brier requires single time") else
                 brier(surv_test, pre_sp = preds_test[, 1], t_star = times),
               "ibs"    = ibs_survmat(surv_test, preds_test, times),
               "iae"    = iae_survmat(surv_test, preds_test, times),
               "ise"    = ise_survmat(surv_test, preds_test, times),
               stop("Unknown metric: ", metric)
        )
      })) |>
      tidyr::unnest(cols = value)
  })

  # Fit final survmetalearner on full data
  tf <- terms(formula, data = data)
  outcome <- attr(tf, "variables")[[2]]
  time_col <- as.character(outcome[[2]])
  status_expr <- outcome[[3]]

  if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
    status_col <- as.character(status_expr[[2]])
    event_value <- eval(status_expr[[3]], data)
    status_vector <- as.integer(data[[status_col]] == event_value)
  } else {
    status_col <- as.character(status_expr)
    event_value <- 1
    status_vector <- data[[status_col]]
  }

  meta_model <- fit_survmetalearner(
    base_preds = base_preds,
    time = data[[time_col]],
    status = status_vector,
    times = times,
    base_models = base_models,
    formula = formula,
    data = data
  )

  structure(list(
    model = meta_model,
    cv_results = results,
    summary = results |> dplyr::group_by(metric) |> dplyr::summarise(mean = mean(value), sd = sd(value)),
    folds = folds,
    metrics = metrics
  ), class = "cv_survmetalearner_result")
}




cv_result <- cv_survmetalearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(200, 500, 800),
  base_models = base_models,
  base_preds = base_preds,
  folds = 3,
  metrics = c("cindex", "ibs", "iae")
)


summary(cv_result$model)
cv_result$summary
plot_survmetalearner_weights(cv_result$model)
score_survmodel(cv_result$model, times = c(200, 500, 800), metrics = c("cindex", "ibs"))


