# ---- cv_survlearner ----
cv_survlearner <- function(formula, data,
                           fit_fun, pred_fun,
                           times,
                           metrics = c("cindex", "ibs"),
                           folds = 5,
                           seed = 123,
                           verbose = TRUE,
                           ...) {

  # warn if '.' is used in formul
  if ("." %in% all.vars(update(formula, . ~ 0))) {
    warning("Please avoid using '.' in the formula. Specify all predictors explicitly.")
  }

  # parse formula
  tf <- terms(formula, data = data)
  outcome <- attr(tf, "variables")[[2]]

  ## need to control
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
    event_value <- 1
    recode_status <- FALSE
  }

  # drop NA in formula variables
  all_vars <- all.vars(formula)
  n_before <- nrow(data)
  data <- tidyr::drop_na(data, dplyr::all_of(all_vars))
  n_after <- nrow(data)


  ## not sure to keep it here, it is informative but it depend on the workload
  if (verbose && n_after < n_before) {
    message("Dropped ", n_before - n_after, " rows with missing values (from ", n_before, " to ", n_after, ").")
  }

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

      # extract outcome for evaluation
      tf <- terms(formula, data = test)
      outcome <- attr(tf, "variables")[[2]]
      time_col <- as.character(outcome[[2]])
      status_expr <- outcome[[3]]

      if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
        status_col <- as.character(status_expr[[2]])
        event_value <- eval(status_expr[[3]], test)
        status_vector <- as.integer(test[[status_col]] == event_value)
      } else {
        status_col <- as.character(status_expr)
        event_value <- 1
        status_vector <- test[[status_col]]
      }

      surv_obj <- survival::Surv(time = test[[time_col]], event = status_vector)

      tibble::tibble(metric = metrics) |>
        dplyr::mutate(value = purrr::map(metric, function(metric) {
          switch(metric,
                 "cindex" = Cindex_survmat(surv_obj, predicted = pred, t_star = max(times)),
                 ## keep it for the moment but likely it should be removed
                 "brier"  = {
                   if (length(times) != 1) stop("Brier requires a single time point.")
                   Brier(surv_obj, pre_sp = pred[, 1], t_star = times)
                 },
                 "ibs"    = Brier_IBS_survmat(surv_obj, sp_matrix = pred, times = times),
                 "iae"    = IAEISE_survmat(surv_obj, sp_matrix = pred, times = times)["IAE"],
                 "ise"    = IAEISE_survmat(surv_obj, sp_matrix = pred, times = times)["ISE"],
                 stop("Unknown metric: ", metric)
          )
        })) |>
        tidyr::unnest(cols = value) |>
        dplyr::mutate(splits = list(split), id = id, fold = fold) |>
        dplyr::relocate(splits, id, fold)
    }
  )
}


library(survival)
library(timereg)
library(rsample)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# define the  model functions if not already
# fit_aareg <- ...
# predict_aareg <- ...

# Example CV evaluation
cv_results <- cv_survlearner(
  formula = Surv(time, status) ~ karno,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs", "iae", "ise"),
  folds = 50
)

# View output
print(cv_results)

# Summarize performance
cv_summary <- function(cv_results) {
  cv_results |>
    group_by(metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      n    = dplyr::n(),
      se   = sd/sqrt(n),
      lower = (mean - 1.96 * se),
      upper = (mean + 1.96 * se),
      .groups = "drop"
    )
}

# Plotting
cv_plot <- function(cv_results) {
  ggplot(cv_results, aes(x = metric, y = value)) +
    geom_boxplot(fill = "skyblue", outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    labs(title = "Cross-Validation Performance", x = "Metric", y = "Value") +
    theme_minimal()
}

# Results
cv_summary(cv_results)
cv_plot(cv_results)
