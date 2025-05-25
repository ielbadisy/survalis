fit_xgboost <- function(formula, data,
                        booster = "gbtree",
                        objective = "survival:aft",
                        aft_loss_distribution = "extreme",
                        aft_loss_distribution_scale = 1,
                        nrounds = 100,
                        ...) {
  stopifnot(requireNamespace("xgboost", quietly = TRUE))

  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(formula, data = mf)[, -1, drop = FALSE]  # remove intercept

  y_time <- y[, "time"]
  y_event <- y[, "status"] == 1

  dmat <- xgboost::xgb.DMatrix(data = x)
  if (objective == "survival:aft") {
    y_lower <- y_time
    y_upper <- y_time
    y_upper[!y_event] <- Inf
    xgboost::setinfo(dmat, "label_lower_bound", y_lower)
    xgboost::setinfo(dmat, "label_upper_bound", y_upper)
  } else if (objective == "survival:cox") {
    xgboost::setinfo(dmat, "label", y_time)
  } else {
    stop("Unsupported objective.")
  }

  params <- list(
    booster = booster,
    objective = objective,
    aft_loss_distribution = aft_loss_distribution,
    aft_loss_distribution_scale = aft_loss_distribution_scale,
    eval_metric = "cox-nloglik"
  )

  model <- xgboost::xgboost(
    params = params,
    data = dmat,
    nrounds = nrounds,
    verbose = 0,
    ...
  )

  structure(list(
    model = model,
    learner = "xgboost",
    formula = formula,
    data = data,
    time = all.vars(formula)[1],
    status = all.vars(formula)[2],
    objective = objective,
    dist = aft_loss_distribution,
    scale = aft_loss_distribution_scale,
    terms = terms(formula)
  ), class = "mlsurv_model", engine = "xgboost")
}


predict_xgboost <- function(object, newdata, times = NULL) {
  stopifnot(object$learner == "xgboost")
  stopifnot(requireNamespace("xgboost", quietly = TRUE))

  mf <- model.frame(object$terms, newdata)
  x_new <- model.matrix(object$terms, data = mf)[, -1, drop = FALSE]

  preds <- predict(object$model, newdata = x_new, outputmargin = TRUE)

  if (is.null(times)) {
    return(preds)  # return linear predictor
  }

  log_times <- outer(rep(1, length(preds)), log(times))
  quants <- (log_times - preds) / object$scale

  surv_probs <- switch(object$dist,
                       "normal" = 1 - pnorm(quants),
                       "extreme" = exp(-exp(quants)),
                       "logistic" = 1 / (1 + exp(quants)),
                       stop("Unsupported distribution.")
  )

  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(x_new)))
  colnames(surv_probs) <- paste0("t=", times)

  as.data.frame(surv_probs)
}



mod_xgb <- fit_xgboost(Surv(time, status) ~ age + karno + celltype, data = veteran)

predict_xgboost(mod_xgb, newdata = veteran[1:5, ], times = c(100, 200, 300))
