# ---- cv-survmodel.R ----
cv_survmodel <- function(data, time_col, status_col, fit_fun, pred_fun, times,
  metrics = c("cindex", "ibs"), folds = 5, event_value = 1, seed = 123, ...) {
set.seed(seed)
rs <- rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col))

results <- purrr::map_dfr(seq_along(rs$splits), function(i) {
split <- rs$splits[[i]]
train_data <- rsample::analysis(split)
test_data  <- rsample::assessment(split)

# Fit model and predict
# Build explicit formula: Surv(time, status == event) ~ x1 + x2 + ...
vars <- setdiff(colnames(train_data), c(time_col, status_col))
rhs <- paste(vars, collapse = " + ")
lhs <- paste0("Surv(", time_col, ", ", status_col, " == ", event_value, ")")
fml <- as.formula(paste(lhs, "~", rhs))

# Fit model
model <- fit_fun(formula = fml, data = train_data, ...)
pred  <- pred_fun(model, newdata = test_data, times = times)

# Evaluate metrics
metric_vals <- evaluate_survmodel(
data = test_data,
time_col = time_col,
status_col = status_col,
event_value = event_value,
sp_matrix = pred,
times = times,
metrics = metrics
)

metric_vals |>
dplyr::mutate(fold = i)
})

return(results)
}


cv_results <- cv_survmodel(
  data = sTRACE,
  time_col = "time",
  status_col = "status",
  event_value = 9,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(1, 3, 5),
  metrics = c("cindex", "ibs"),
  folds = 5
)

print(cv_results)



########



cv_survmodel <- function(data, time_col, status_col,
  fit_fun, pred_fun, times,
  metrics = c("cindex", "ibs"),
  folds = 5, event_value = 1,
  seed = 123, ...) {
set.seed(seed)

# Create v-folds
rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col)) |>
dplyr::mutate(
fold = dplyr::row_number(),
results = purrr::map(splits, function(split) {
train_data <- rsample::analysis(split)
test_data  <- rsample::assessment(split)

# Explicit formula with variable names
predictors <- setdiff(colnames(train_data), c(time_col, status_col))
formula_str <- paste0("Surv(", time_col, ", ", status_col, " == ", event_value, ") ~ ",
       paste(predictors, collapse = " + "))
fml <- as.formula(formula_str)

# Fit and predict
model <- fit_fun(formula = fml, data = train_data, ...)
pred  <- pred_fun(model, newdata = test_data, times = times)

# Evaluate metrics
evaluate_survmodel(
data = test_data,
time_col = time_col,
status_col = status_col,
event_value = event_value,
sp_matrix = pred,
times = times,
metrics = metrics
)
})
) |>
tidyr::unnest(results)
}



cv_results <- cv_survmodel(
  data = sTRACE,
  time_col = "time",
  status_col = "status",
  event_value = 9,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(1, 3, 5),
  metrics = c("cindex", "ibs"),
  folds = 5
)

print(cv_results)
cv_results <- cv_survmodel(
  data = sTRACE,
  time_col = "time",
  status_col = "status",
  event_value = 9,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(1, 3, 5),
  metrics = c("cindex", "ibs"),
  folds = 5
)

print(cv_results)






#########



cv_summary <- cv_results |>
  dplyr::group_by(metric) |>
  dplyr::summarise(
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

print(cv_summary)




### 


cv_summary <- function(cv_results) {
  cv_results |>
    dplyr::group_by(metric) |>
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      n    = dplyr::n(),
      se   = sd / sqrt(n),
      lower = mean - 1.96 * se,
      upper = mean + 1.96 * se,
      .groups = "drop"
    )
}


cv_summary(cv_results)




cv_plot <- function(cv_results) {
  ggplot2::ggplot(cv_results, ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.4, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.6) +
    ggplot2::labs(
      title = "Cross-Validation Performance by Metric",
      y = "Metric Value", x = "Metric"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}



cv_plot(cv_results)



cv_plot_faceted <- function(cv_results) {
  if (!"learner" %in% colnames(cv_results)) {
    stop("cv_results must include a 'learner' column to facet by.")
  }

  ggplot2::ggplot(cv_results, ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_boxplot(fill = "skyblue", alpha = 0.3, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.6) +
    ggplot2::facet_wrap(~ learner) +
    ggplot2::labs(
      title = "Cross-Validation Performance by Learner and Metric",
      x = "Metric", y = "Value"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}



learners <- list(
  aareg1 = list(fit = fit_aareg, pred = predict_aareg),
  aareg12   = list(fit = fit_aareg, pred = predict_aareg)
)

multi_cv <- purrr::imap_dfr(learners, function(learner_fns, name) {
  cv_survmodel(
    data = sTRACE,
    time_col = "time",
    status_col = "status",
    event_value = 9,
    fit_fun = learner_fns$fit,
    pred_fun = learner_fns$pred,
    times = c(1, 3, 5),
    metrics = c("cindex", "ibs"),
    folds = 5
  ) |> dplyr::mutate(learner = name)
})



cv_plot_faceted(multi_cv)
