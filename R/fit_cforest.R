#' Fit a Conditional Inference Survival Forest
#'
#' Fits a survival forest model using the `party::cforest()` implementation
#' of conditional inference trees for right-censored survival data.
#'
#' @param formula A `Surv()` survival formula specifying the outcome and predictors.
#' @param data A data frame containing the variables in the model.
#' @param teststat Character string specifying the test statistic to use 
#'   (default = `"quad"`). See [party::cforest_control()].
#' @param testtype Character string specifying the type of test 
#'   (default = `"Univariate"`).
#' @param mincriterion Numeric, the value of the test statistic that must be
#'   exceeded for a split to be performed (default = `0`).
#' @param ntree Integer, number of trees to grow (default = `500`).
#' @param mtry Integer, number of variables randomly selected at each node
#'   (default = `5`).
#' @param replace Logical, whether sampling of cases is with replacement
#'   (default = `TRUE`).
#' @param fraction Proportion of samples to draw if `replace = FALSE`
#'   (default = `0.632`).
#' @param ... Additional arguments passed to [party::cforest_control()].
#'
#' @return An object of class `"mlsurv_model"`, containing:
#'   \item{model}{The fitted `cforest` model.}
#'   \item{learner}{Character string `"cforest"`.}
#'   \item{formula}{The model formula.}
#'   \item{data}{The training data.}
#'
#' @details This function wraps `party::cforest()` to fit a conditional inference
#'   forest for survival analysis, returning an `mlsurv_model` object compatible
#'   with `predict_cforest()` and the `survalis` framework.
#'
#' @examples
#' mod <- fit_cforest(Surv(time, status) ~ age + celltype + karno, data = veteran)
#' head(predict_cforest(mod, newdata = veteran[1:5, ], times = c(100, 200)))
#'
#' @export

fit_cforest <- function(formula, data,
                        teststat = "quad",
                        testtype = "Univariate",
                        mincriterion = 0,
                        ntree = 500,
                        mtry = 5,
                        replace = TRUE,
                        fraction = 0.632,
                        ...) {
  stopifnot(requireNamespace("party", quietly = TRUE))

  model <- suppressWarnings(
    party::cforest(
      formula = formula,
      data = data,
      controls = party::cforest_control(
        teststat = teststat,
        testtype = testtype,
        mincriterion = mincriterion,
        ntree = ntree,
        mtry = mtry,
        replace = replace,
        fraction = fraction,
        ...
      )
    )
  )

  structure(
    list(
      model = model,
      learner = "cforest",
      formula = formula,
      data = data
    ),
    class = "mlsurv_model",
    engine = "party"
  )
}


#' Predict Survival Probabilities from a Conditional Inference Survival Forest
#'
#' Generates predicted survival probabilities at specified time points from a
#' fitted `cforest` survival model.
#'
#' @param object An `mlsurv_model` object returned by [fit_cforest()].
#' @param newdata A data frame containing the predictor variables for prediction.
#' @param times A numeric vector of time points at which to estimate survival
#'   probabilities.
#' @param ... Not used, included for compatibility.
#'
#' @return A data frame of survival probabilities with one row per observation
#'   in `newdata` and one column per requested time point (named `"t=<time>"`).
#'
#' @details Survival curves are extracted from the fitted `cforest` model and
#'   linearly interpolated to the requested time points.
#'
#' @examples
#' mod <- fit_cforest(Surv(time, status) ~ age + celltype + karno, data = veteran)
#' predict_cforest(mod, newdata = veteran[1:5, ], times = c(100, 200, 300))
#'
#' @export

predict_cforest <- function(object, newdata, times, ...) {



  if (!is.null(object$learner) && object$learner != "cforest") {
    warning("Object passed to predict_cforest() may not come from fit_cforest().")
  }



  stopifnot(requireNamespace("party", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  # get survival predictions from cforest
  survlist <- suppressMessages(suppressWarnings(
    predict(object$model, newdata = newdata, type = "prob")
  ))

  # interpolate at requested times
  survmat <- t(sapply(survlist, function(sfit) {
    s_time <- sfit$time
    s_surv <- sfit$surv

    if (is.null(s_time) || is.null(s_surv) || length(s_surv) == 0) {
      return(rep(NA, length(times)))
    }

    approx(x = s_time, y = s_surv, xout = times, method = "linear", rule = 2)$y
  }))

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}

#' Tune a Conditional Inference Survival Forest
#'
#' Performs cross-validated hyperparameter tuning for a conditional inference
#' survival forest using the `fit_cforest()` and `predict_cforest()` functions.
#'
#' @param formula A `Surv()` survival formula specifying the outcome and predictors.
#' @param data A data frame containing the variables in the model.
#' @param times A numeric vector of time points at which to evaluate performance.
#' @param param_grid A data frame or list specifying the grid of hyperparameter
#'   values to evaluate. Columns should include `ntree`, `mtry`, `mincriterion`,
#'   and `fraction`.
#' @param metrics A character vector of performance metrics to compute 
#'   (default = `c("cindex", "ibs")`).
#' @param folds Integer, number of cross-validation folds (default = `5`).
#' @param seed Integer, random seed for reproducibility.
#' @param refit_best Logical, if `TRUE` refits a model with the best-performing
#'   hyperparameters and returns it; otherwise returns a summary table of results.
#' @param ... Additional arguments passed to [fit_cforest()].
#'
#' @return If `refit_best = FALSE`, a data frame of mean cross-validation scores
#'   for each hyperparameter combination (class `"tuned_surv"`).  
#'   If `refit_best = TRUE`, an `mlsurv_model` object fitted with the optimal
#'   hyperparameters.
#'
#' @details Cross-validation is performed using `cv_survlearner()` and results
#'   are sorted in descending order of the first metric specified in `metrics`.
#'
#' @examples
#' grid <- expand.grid(
#'   ntree = c(50, 100),
#'   mtry = c(2, 4),
#'   mincriterion = c(0, 0.95),
#'   fraction = c(0.632)
#' )
#' res <- tune_cforest(Surv(time, status) ~ age + celltype + karno,
#'   data = veteran,
#'   times = c(100, 200),
#'   param_grid = grid,
#'   folds = 3
#' )
#' print(res)
#'
#' @export

tune_cforest <- function(formula, data, times,
                         param_grid = expand.grid(
                           ntree = c(100, 300),
                           mtry = c(2, 3),
                           mincriterion = c(0, 0.95),
                           fraction = c(0.5, 0.632)
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5, seed = 123,
                         refit_best = FALSE, ...) {

  results <- purrr::pmap_dfr(param_grid, function(ntree, mtry, mincriterion, fraction) {
    cv_results <- cv_survlearner(
      formula = formula,
      data = data,
      fit_fun = fit_cforest,
      pred_fun = predict_cforest,
      times = times,
      metrics = metrics,
      folds = folds,
      seed = seed,
      ntree = ntree,
      mtry = mtry,
      mincriterion = mincriterion,
      fraction = fraction,
      ...
    )

    summary <- cv_summary(cv_results)

    tibble::tibble(ntree = ntree,
                   mtry = mtry,
                   mincriterion = mincriterion,
                   fraction = fraction) |>
      dplyr::bind_cols(
        tidyr::pivot_wider(summary[, c("metric", "mean")],
                           names_from = metric,
                           values_from = mean)
      )
  })

  if (metrics[1] %in% c("cindex", "auc", "accuracy")) {
    results <- results |> dplyr::arrange(dplyr::desc(!!rlang::sym(metrics[1])))
  } else {
    results <- results |> dplyr::arrange(!!rlang::sym(metrics[1]))
  }

  if (!refit_best) {
    class(results) <- c("tuned_surv", class(results))
    attr(results, "metrics") <- metrics
    attr(results, "formula") <- formula
    return(results)
  } else {
    best_row <- results[1, ]
    best_model <- fit_cforest(
      formula = formula,
      data = data,
      ntree = best_row$ntree,
      mtry = best_row$mtry,
      mincriterion = best_row$mincriterion,
      fraction = best_row$fraction,
      ...
    )
    return(best_model)
  }
}
