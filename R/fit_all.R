# Lists all learners with matching `fit_*` and `predict_*` functions
survlearners <- function() {
  fits <- ls(pattern = "^fit_", envir = .GlobalEnv)
  learners <- sub("^fit_", "", fits)
  learners[paste0("predict_", learners) %in% ls(pattern = "^predict_", envir = .GlobalEnv)]
}



# define a tibble of survival learners (specific to benchmarking)
define_survlearners <- function(...) {
  learner_names <- c(...)
  if (length(learner_names) == 0) {
    learner_names <- survlearners()
  }
  tibble::tibble(
    learner_name = learner_names,
    fit_fun = purrr::map(learner_names, ~ get(paste0("fit_", .x), mode = "function")),
    pred_fun = purrr::map(learner_names, ~ get(paste0("predict_", .x), mode = "function"))
  )
}


# benchmark multiple survival learners using cv
benchmark_survlearners <- function(learners,
                                   formula,
                                   data,
                                   times,
                                   metrics = c("cindex", "ibs"),
                                   folds = 5,
                                   seed = 123,
                                   ...) {
  learners %>%
    dplyr::mutate(results = purrr::map2(fit_fun, pred_fun, function(fit, pred) {
      cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit,
        pred_fun = pred,
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        ...
      )
    })) %>%
    dplyr::select(learner_name, results) %>%
    tidyr::unnest(results)
}



## TEST

## first we need to define the learners, this can be done manualy through building and tibble or using helpers

# define models and prediction functions manually
learners <- tibble::tibble(
  learner_name = c("aareg", "aftgee"),
  fit_fun = list(fit_aareg, fit_aftgee),
  pred_fun = list(predict_aareg, predict_aftgee)
)


# define models and prediction functions through helpers

learners <- define_survlearners("aareg", "aftgee")

benchmark_results <- benchmark_survlearners(
  learners = learners,
  formula = Surv(Time, status) ~ x1 + x2,
  data = dat,
  times = c(10, 20, 40),
  metrics = c("cindex", "ibs"),
  folds = 5
)

print(benchmark_results)




#****** summarise benchmark results across folds

summarise_benchmark <- function(benchmark_results) {
  benchmark_results |>
    dplyr::group_by(learner_name, metric) |>
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


# plot benchmark results

plot_benchmark <- function(benchmark_results) {
  ggplot2::ggplot(benchmark_results, ggplot2::aes(x = learner_name, y = value)) +
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


## TEST
summarise_benchmark(benchmark_results)
plot_benchmark(benchmark_results)

#***************** helper to select the best survleanrer

# select the best survival learner per metric

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
    dplyr::group_by(learner_name, metric) |>
    dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

  best <- summary |>
    dplyr::filter(if (maximize) value == max(value) else value == min(value))

  return(best)
}

### TO POLISH table names (maybe unify columns names like learner_name -> survlearner)
best_survlearner(benchmark_results, "ibs")

#************** for more

learners <- tibble::tibble(
  learner_name = c("aareg", "aftgee", "mboost"),
  fit_fun = list(fit_aareg, fit_aftgee, fit_mboost),
  pred_fun = list(predict_aareg, predict_aftgee, predict_mboost)
)

benchmark_results <- benchmark_survlearners(
  learners = learners,
  formula = Surv(Time, status) ~ x1 + x2,
  data = dat,
  times = c(10, 20, 40),
  metrics = c("cindex", "ibs"),
  folds = 5
)

