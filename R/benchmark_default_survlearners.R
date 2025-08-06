
benchmark_default_survlearners <- function(formula, data, learners, times,
                                   metrics = c("cindex", "ibs"),
                                   folds = 5, seed = 123,
                                   verbose = FALSE,
                                   suppress_errors = TRUE,
                                   ...) {
  stopifnot(is.character(learners), length(learners) > 0)

  results_list <- lapply(learners, function(learner) {
    fit_fun_name  <- paste0("fit_", learner)
    pred_fun_name <- paste0("predict_", learner)

    if (!exists(fit_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", fit_fun_name))
      return(NULL)
    }
    if (!exists(pred_fun_name, mode = "function")) {
      warning(sprintf("Function '%s' not found.", pred_fun_name))
      return(NULL)
    }

    fit_fun  <- get(fit_fun_name)
    pred_fun <- get(pred_fun_name)

    if (verbose) message("Benchmarking learner: ", learner)

    result <- try(
      cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_fun,
        pred_fun = pred_fun,
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        verbose = verbose,
        ...
      ),
      silent = suppress_errors
    )

    if (inherits(result, "try-error")) {
      if (!suppress_errors) stop(result)
      warning(sprintf("Learner '%s' failed: %s", learner, conditionMessage(attr(result, "condition"))))
      return(NULL)
    }

    result$learner <- learner
    result
  })

  # Combine successful results
  results_combined <- dplyr::bind_rows(results_list)
  if (nrow(results_combined) == 0) {
    stop("All learners failed or returned empty results.")
  }

  results_combined
}



learners <- c("aalen", "aftgee", "blackboost", "ranger", "orsf", "glmnet", "flexsurvreg", "bnnsurv", "coxaalen", "stpm2", "rpart", "rsf", "survsvm", "xgboost", "coxph", "survdnn")


results <- benchmark_default_survlearners(
  formula = Surv(time, status) ~ age + trt + karno,
  data = veteran,
  learners = learners,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs", "iae", "ise"),
  verbose = TRUE
)

dplyr::glimpse(results)


summarise_benchmark <- function(benchmark_results) {
  benchmark_results |>
    dplyr::group_by(learner, metric) |>
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


summarise_benchmark(results)


plot_benchmark <- function(benchmark_results) {
  ggplot2::ggplot(benchmark_results, ggplot2::aes(x = learner, y = value)) +
    ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.4, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
    ggplot2::facet_wrap(~ metric, scales = "free_y") +
    ggplot2::labs(
      title = "Cross-Validated Performance of Survival Learners",
      x = "Learner",
      y = "Metric Value"
    ) +
    ggplot2::theme_minimal(base_size = 14)
}


plot_benchmark(results)


summarize_benchmark_results <- function(results, digits = 3) {

  stopifnot(is.data.frame(results), all(c("learner", "metric", "value") %in% colnames(results)))

  results %>%
    dplyr::group_by(learner, metric) %>%
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      summary = sprintf("%.*f ± %.*f", digits, mean, digits, sd)
    ) %>%
    tidyr::pivot_wider(
      id_cols = learner,
      names_from = metric,
      values_from = summary
    ) %>%
    dplyr::arrange(learner)
}

summary_tbl <- summarize_benchmark_results(results)
print(summary_tbl)





best_survlearner <- function(benchmark_results, metric, maximize = NULL) {
  if (!metric %in% benchmark_results$metric) {
    stop("Metric not found in results: ", metric)
  }

  if (is.null(maximize)) {
    # default: maximize cindex, minimize others
    maximize <- !(metric %in% c("ibs", "brier", "iae", "ise"))
  }

  summary <- benchmark_results |>
    dplyr::filter(metric == !!metric) |>
    dplyr::group_by(learner, metric) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  best <- summary |>
    dplyr::filter(if (maximize) value == max(value) else value == min(value))

  return(best)
  }

### TO POLISH table names (maybe unify columns names like learner_name -> survlearner)
best_survlearner(results, "cindex")
