
#' Benchmark Multiple Survival Learners (Cross-Validation Wrapper)
#'
#' Runs \code{cv_survlearner()} for a set of learner names (e.g., \code{"ranger"},
#' \code{"coxph"}) by dynamically dispatching \code{fit_<learner>} and
#' \code{predict_<learner>} functions. Returns the row‑bound CV results across
#' all requested learners.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ x1 + x2 + ...}.
#' @param data A data frame containing the variables in \code{formula}.
#' @param learners Character vector of learner ids (without prefixes), e.g.
#'   \code{c("ranger","coxph","glmnet")}. For each id, \code{fit_<id>} and
#'   \code{predict_<id>} must exist.
#' @param times Numeric vector of evaluation time points for survival predictions.
#' @param metrics Character vector of metrics to compute in CV (e.g., \code{"cindex"},
#'   \code{"ibs"}, \code{"iae"}, \code{"ise"}). The metric semantics are those of
#'   \code{cv_survlearner()}.
#' @param folds Integer number of CV folds (default \code{5}).
#' @param seed Integer random seed for fold generation.
#' @param verbose Logical; if \code{TRUE}, prints progress per learner.
#' @param suppress_errors Logical; if \code{TRUE} (default) errors from individual
#'   learners are caught and reported via \code{warning()} and benchmarking continues.
#'   If \code{FALSE}, the first error is thrown.
#' @param ... Additional arguments forwarded to each learner's \code{fit_*()}.
#'
#' @details
#' Learners are run independently using identical CV splits and scoring settings.
#' Any learner whose \code{fit_*()} or \code{predict_*()} function is missing will
#' be skipped with a warning. At least one learner must complete successfully or an
#' error is raised.
#'
#' @return A data frame of CV results (as returned by \code{cv_survlearner()})
#' with an extra column \code{learner} identifying the source learner.
#'
#' @examples
#' # learners <- c("coxph","ranger","glmnet")
#' # res <- benchmark_default_survlearners(
#' #   Surv(time, status) ~ age + karno + trt, veteran,
#' #   learners = learners, times = c(100,200,300),
#' #   metrics = c("cindex","ibs"), folds = 3, seed = 42
#' # )
#' # dplyr::glimpse(res)
#'
#' @seealso [cv_survlearner()], [plot_benchmark()], [summarise_benchmark()]
#' @export

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



#' Summarise Benchmark Results (Mean  SD with Wald CI)
#'
#' Aggregates cross‑validated benchmark results by learner and metric, reporting
#' mean, standard deviation, standard error, and an approximate 95% Wald
#' confidence interval.
#'
#' @param benchmark_results A data frame produced by
#'   \code{benchmark_default_survlearners()}, containing at least
#'   \code{learner}, \code{metric}, and \code{value}.
#'
#' @return A tibble with columns \code{learner}, \code{metric}, \code{mean},
#' \code{sd}, \code{n}, \code{se}, \code{lower}, \code{upper}.
#'
#' @examples
#' # summ <- summarise_benchmark(res)
#' # summ
#'
#' @seealso [benchmark_default_survlearners()], [plot_benchmark()],
#'   [summarize_benchmark_results()]
#' @export

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



#' Plot Benchmark Distributions Across Learners
#'
#' Produces box‑and‑jitter plots of CV metric values per learner, faceted by
#' metric for quick visual comparison.
#'
#' @param benchmark_results A data frame from
#'   \code{benchmark_default_survlearners()} with columns \code{learner},
#'   \code{metric}, and \code{value}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' # plot_benchmark(res)
#'
#' @seealso [benchmark_default_survlearners()], [summarise_benchmark()]
#' @export

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



#' Compact Table of Mean  SD by Learner and Metric
#'
#' Creates a wide table summarizing each learner's performance as
#' \code{mean  sd} per metric, suitable for reporting.
#'
#' @param results A data frame with columns \code{learner}, \code{metric},
#'   and \code{value} (e.g., output of \code{benchmark_default_survlearners()}).
#' @param digits Integer number of decimal places in the formatted summary
#'   (default \code{3}).
#'
#' @return A wide tibble with one row per learner and one column per metric,
#'   containing formatted strings \code{"mean  sd"}.
#'
#' @examples
#' # tbl <- summarize_benchmark_results(res, digits = 2)
#' # tbl
#'
#' @seealso [summarise_benchmark()], [benchmark_default_survlearners()]
#' @export

summarize_benchmark_results <- function(results, digits = 3) {

  stopifnot(is.data.frame(results), all(c("learner", "metric", "value") %in% colnames(results)))

  results |>
    dplyr::group_by(learner, metric) |>
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      summary = sprintf("%.*f  %.*f", digits, mean, digits, sd)
    ) |>
    tidyr::pivot_wider(
      id_cols = learner,
      names_from = metric,
      values_from = summary
    ) |>
    dplyr::arrange(learner)
}





#' Select the Best Survival Learner by a Given Metric
#'
#' Extracts the top‑performing learner(s) under a chosen metric from benchmark
#' results, using the average value across folds.
#'
#' @param benchmark_results A data frame with columns \code{learner},
#'   \code{metric}, and \code{value}.
#' @param metric Character name of the metric to optimize (e.g., \code{"cindex"},
#'   \code{"ibs"}). Must exist in \code{benchmark_results$metric}.
#' @param maximize Logical; whether to maximize the metric. If \code{NULL}
#'   (default), the function maximizes for concordance‑like metrics
#'   (\code{"cindex"}) and minimizes for error‑like metrics
#'   (\code{"ibs"}, \code{"brier"}, \code{"iae"}, \code{"ise"}).
#'
#' @return A tibble with columns \code{learner}, \code{metric}, and the selected
#'   average \code{value} for the best learner(s). Ties are returned as multiple rows.
#'
#' @examples
#' # best_survlearner(res, metric = "cindex")
#' # best_survlearner(res, metric = "ibs")  # minimized by default
#'
#' @seealso [benchmark_default_survlearners()], [summarise_benchmark()]
#' @export

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


