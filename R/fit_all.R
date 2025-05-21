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



