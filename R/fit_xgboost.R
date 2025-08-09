
fit_xgboost <- function(formula, data,
                        objective = "survival:aft",
                        nrounds = 100,
                        max_depth = 3,
                        eta = 0.1,
                        aft_loss_distribution = "normal",
                        aft_loss_distribution_scale = 1.0,
                        ...) {

  stopifnot(requireNamespace("xgboost", quietly = TRUE))


  mf <- model.frame(formula, data)
  y <- model.response(mf)

  time <- y[, 1]
  status <- y[, 2]
  x <- model.matrix(attr(mf, "terms"), mf)[, -1, drop = FALSE]

  if (objective == "survival:aft") {
    y_lower <- ifelse(status == 1, time, -Inf)
    y_upper <- ifelse(status == 1, time, Inf)
    ymat <- cbind(y_lower, y_upper)

    dtrain <- xgboost::xgb.DMatrix(data = x, label_lower_bound = ymat[, 1],
                                   label_upper_bound = ymat[, 2])
    params <- list(
      objective = "survival:aft",
      max_depth = max_depth,
      eta = eta,
      aft_loss_distribution = aft_loss_distribution,
      aft_loss_distribution_scale = aft_loss_distribution_scale
    )
  } else if (objective == "survival:cox") {
    dtrain <- xgboost::xgb.DMatrix(data = x, label = time)
    dtrain$setinfo("label", time)
    dtrain$setinfo("label_lower_bound", time)
    dtrain$setinfo("label_upper_bound", time)
    dtrain$setinfo("label", time)
    dtrain$setinfo("event", status)

    params <- list(
      objective = "survival:cox",
      max_depth = max_depth,
      eta = eta
    )
  } else {
    stop("Unsupported objective")
  }

  bst <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    verbose = 0
  )

  structure(list(
    model = bst,
    data = data,
    learner = "xgboost",
    formula = formula,
    xnames = colnames(x),
    objective = objective
  ), class = "mlsurv_model", engine = "xgboost")
}


predict_xgboost <- function(object, newdata, times) {


  if (!is.null(object$learner) && object$learner != "xgboost") {warning("Object passed to predict_xgboost() may not come from fit_xgboost().")}

  stopifnot(requireNamespace("xgboost", quietly = TRUE))
  

  xmat <- model.matrix(terms(object$formula), newdata)[, object$xnames, drop = FALSE]

  if (object$objective == "survival:aft") {
    mu <- predict(object$model, newdata = xmat)
    sigma <- object$model$params$aft_loss_distribution_scale %||% 1.0
    dist <- object$model$params$aft_loss_distribution %||% "normal"

    survmat <- switch(dist,
                         "normal" = outer(mu, times, function(m, t) pnorm((log(t) - m) / sigma, lower.tail = FALSE)),
                         "logistic" = outer(mu, times, function(m, t) plogis((log(t) - m) / sigma, lower.tail = FALSE)),
                         "extreme" = outer(mu, times, function(m, t) exp(-exp((log(t) - m) / sigma))),
                         stop("Unsupported AFT distribution")
    )
  } else if (object$objective == "survival:cox") {
    lp <- predict(object$model, newdata = xmat)
    basehaz <- survfit(coxph(Surv(rep(1, length(lp)), rep(1, length(lp))) ~ lp))
    bh_fun <- stepfun(basehaz$time, c(0, basehaz$cumhaz))
    survmat <- sapply(times, function(t) exp(-bh_fun(t) * exp(lp)))
    survmat <- t(survmat)
  } else {
    stop("Unsupported objective for prediction")
  }

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}

tune_xgboost <- function(formula, data, times,
                         param_grid,
                         metrics = c("cindex", "ibs"),
                         folds = 5,
                         seed = 42,
                         refit_best = TRUE) {
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(tibble)

  # Generate grid of hyperparameters
  grid_df <- tidyr::crossing(!!!param_grid)

  # Perform CV for each grid combination
  results <- purrr::pmap_dfr(grid_df, function(...) {
    params <- list(...)

    fit_fun <- function(formula, data) {
      do.call(fit_xgboost, c(list(formula = formula, data = data), params))
    }

    res_cv <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_xgboost,
      pred_fun = predict_xgboost,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      verbose = FALSE
    )

    summary_tbl <- res_cv %>%
      group_by(metric) %>%
      summarise(mean = mean(value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = metric, values_from = mean)

    bind_cols(as_tibble(params), summary_tbl)
  })

  # Rank by first metric
  primary_metric <- metrics[1]
  results <- arrange(results, desc(.data[[primary_metric]]))

  if (refit_best && nrow(results) > 0) {
    best_params <- results[1, names(param_grid), drop = FALSE]
    final_model <- do.call(fit_xgboost, c(list(formula = formula, data = data), best_params))
    final_model$data <- data  # required for summary() to work
    return(final_model)
  }

  return(results)
}





mod_xgboost <- fit_xgboost(
  Surv(time, status) ~ age + karno + celltype,
  data = veteran
)

pred_surv <- predict_xgboost(mod_xgboost, newdata = veteran, times = 100)

pred_surv



cv_results_xgb <- cv_survlearner(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  fit_fun = fit_xgboost,
  pred_fun = predict_xgboost,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs", "iae", "ise"),
  folds = 5,
  verbose = TRUE
)



grid <- list(
  nrounds = c(50, 100),
  max_depth = c(3, 6),
  eta = c(0.01, 0.1),
  aft_loss_distribution = c("extreme", "logistic"),
  aft_loss_distribution_scale = c(0.5, 1),
  objective = "survival:aft"
)

tune_res <- tune_xgboost(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 300, 500),
  param_grid = grid,
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

summary(tune_res)



pred_surv <- predict_xgboost(tune_res, newdata = veteran, times = c("cindex", "ibs"))


