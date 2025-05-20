# ---- cv survlearner ----
cv_survlearner <- function(formula, data,
  fit_fun, pred_fun,
  times,
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ...) {
# warn if '.' is used in formula
if ("." %in% all.vars(update(formula, . ~ 0))) {
warning("Please avoid using '.' in the formula. Specify all predictors explicitly.")
}

# parse formula
tf <- terms(formula, data = data)
outcome <- attr(tf, "variables")[[2]]

if (!inherits(outcome, "call") || outcome[[1]] != as.name("Surv")) {
stop("Formula must begin with a Surv(...) outcome.")
}

time_col <- as.character(outcome[[2]])
status_expr <- outcome[[3]]

# flexible status handling
if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
status_col <- as.character(status_expr[[2]])
event_value <- eval(status_expr[[3]], data)
recode_status <- TRUE
} else {
status_col <- as.character(status_expr)
event_value <- 1  # conventionally treat 1 as event
recode_status <- FALSE
}

# drop NA in formula variables
all_vars <- all.vars(formula)
n_before <- nrow(data)
data <- tidyr::drop_na(data, dplyr::all_of(all_vars))
n_after <- nrow(data)

if (n_after < n_before) {
message("Dropped ", n_before - n_after, " rows with missing values (from ", n_before, " to ", n_after, ").")
}

# if status needs recoding -> do it
if (recode_status) {
data[[status_col]] <- as.integer(data[[status_col]] == event_value)
}

set.seed(seed)
folds_data <- rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col))

purrr::pmap_dfr(
list(split = folds_data$splits, id = folds_data$id, fold = seq_len(folds)),
function(split, id, fold) {
train <- rsample::analysis(split)
test  <- rsample::assessment(split)

model <- fit_fun(formula = formula, data = train, ...)
pred  <- pred_fun(model, newdata = test, times = times)

evaluate_survmodel(
data = test,
time_col = time_col,
status_col = status_col,
event_valu# ---- cv survlearner ----
cv_survlearner <- function(formula, data,
  fit_fun, pred_fun,
  times,
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  ...) {
# warn if '.' is used in formula
if ("." %in% all.vars(update(formula, . ~ 0))) {
warning("Please avoid using '.' in the formula. Specify all predictors explicitly.")
}

# parse formula
tf <- terms(formula, data = data)
outcome <- attr(tf, "variables")[[2]]

if (!inherits(outcome, "call") || outcome[[1]] != as.name("Surv")) {
stop("Formula must begin with a Surv(...) outcome.")
}

time_col <- as.character(outcome[[2]])
status_expr <- outcome[[3]]

# flexible status handling
if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
status_col <- as.character(status_expr[[2]])
event_value <- eval(status_expr[[3]], data)
recode_status <- TRUE
} else {
status_col <- as.character(status_expr)
event_value <- 1  # conventionally treat 1 as event
recode_status <- FALSE
}

# drop NA in formula variables
all_vars <- all.vars(formula)
n_before <- nrow(data)
data <- tidyr::drop_na(data, dplyr::all_of(all_vars))
n_after <- nrow(data)

if (n_after < n_before) {
message("Dropped ", n_before - n_after, " rows with missing values (from ", n_before, " to ", n_after, ").")
}

# if status needs recoding -> do it
if (recode_status) {
data[[status_col]] <- as.integer(data[[status_col]] == event_value)
}

set.seed(seed)
folds_data <- rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col))

purrr::pmap_dfr(
list(split = folds_data$splits, id = folds_data$id, fold = seq_len(folds)),
function(split, id, fold) {
train <- rsample::analysis(split)
test  <- rsample::assessment(split)

model <- fit_fun(formula = formula, data = train, ...)
pred  <- pred_fun(model, newdata = test, times = times)

evaluate_survmodel(
data = test,
time_col = time_col,
status_col = status_col,
event_value = 1,
sp_matrix = pred,
times = times,
metrics = metrics
) |>
dplyr::mutate(splits = list(split), id = id, fold = fold) |>
dplyr::relocate(splits, id, fold)
}
)
}


## without specifying the status == value)

cv_survlearner(
  Surv(time, status) ~ karno,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(100, 300, 500),
  folds = 5,
  seed = 123
)


## with specifying the status == value)
cv_survlearner(
  Surv(time, status == 1) ~ karno,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(100),
  folds = 5,
  seed = 123
)


#********** plotting for CV 


## very basic for the moment need more work


## CV summary

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





## plot
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




e = 1,
sp_matrix = pred,
times = times,
metrics = metrics
) |>
dplyr::mutate(splits = list(split), id = id, fold = fold) |>
dplyr::relocate(splits, id, fold)
}
)
}


## without specifying the status == value)

cv_survlearner(
  Surv(time, status) ~ karno,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(100, 300, 500),
  folds = 5,
  seed = 123
)


## with specifying the status == value)
cv_survlearner(
  Surv(time, status == 1) ~ karno,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(100),
  folds = 5,
  seed = 123
)


#********** plotting for CV 


## very basic for the moment need more work


## CV summary

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





## plot
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




