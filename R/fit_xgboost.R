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

colnames(surv_probs) <- paste0("t=", times)

as.data.frame(surv_probs)
}



mod_xgb <- fit_xgboost(Surv(time, status) ~ age + karno + celltype, data = veteran)

predict_xgboost(mod_xgb, newdata = veteran[1:5, ], times = c(100, 200, 300))


#----



tune_xgboost <- function(formula, data, times,
  param_grid = expand.grid(
    booster = c("gbtree"),
    objective = c("survival:aft"),
    aft_loss_distribution = c("extreme", "logistic"),
    aft_loss_distribution_scale = c(0.5, 1, 2),
    nrounds = c(50, 100, 200),
    stringsAsFactors = FALSE
  ),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  refit_best = TRUE,
  ...) {

.higher_is_better <- function(metric) metric %in% c("cindex","auc","accuracy")

results <- purrr::pmap_dfr(param_grid, function(booster, objective,
                           aft_loss_distribution,
                           aft_loss_distribution_scale,
                           nrounds) {
cv_res <- tryCatch({
cv_survlearner(
formula = formula,
data = data,
fit_fun = fit_xgboost,
pred_fun = predict_xgboost,
times = times,
metrics = metrics,
folds = folds,
seed = seed,
booster = booster,
objective = objective,
aft_loss_distribution = aft_loss_distribution,
aft_loss_distribution_scale = aft_loss_distribution_scale,
nrounds = nrounds,
...
)
}, error = function(e) NULL)

if (is.null(cv_res)) {
return(tibble::tibble(
booster, objective, aft_loss_distribution,
aft_loss_distribution_scale, nrounds,
failed = TRUE
))
}

smry <- cv_summary(cv_res)
tibble::tibble(
booster, objective, aft_loss_distribution,
aft_loss_distribution_scale, nrounds,
failed = FALSE
) |>
dplyr::bind_cols(
tidyr::pivot_wider(smry[, c("metric","mean")],
    names_from = metric, values_from = mean)
)
})

primary <- metrics[1]
cand <- results |> dplyr::filter(!failed)
if (nrow(cand) == 0L) stop("All tuning configurations failed.")

ord1 <- if (.higher_is_better(primary)) dplyr::desc(cand[[primary]]) else cand[[primary]]
best_row <- cand[order(ord1), ][1, , drop = FALSE]

if (!refit_best) {
return(results |>
dplyr::arrange(
if (.higher_is_better(primary)) dplyr::desc(.data[[primary]]) else .data[[primary]]
))
}

best_model <- fit_xgboost(
formula = formula,
data = data,
booster = best_row$booster,
objective = best_row$objective,
aft_loss_distribution = best_row$aft_loss_distribution,
aft_loss_distribution_scale = best_row$aft_loss_distribution_scale,
nrounds = best_row$nrounds,
...
)

attr(best_model, "tuning_results") <- results
attr(best_model, "best_params") <- as.list(best_row)
best_model
}




mod_xgb_best <- tune_xgboost(
  formula = Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 3,
  refit_best = TRUE
)

# Works directly with predict_*
predict_xgboost(mod_xgb_best, newdata = veteran[1:5, ], times = c(100, 200, 300))



summary(mod_xgb_best)
