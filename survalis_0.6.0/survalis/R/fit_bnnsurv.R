

#' Fit a Survival SVM Model (survivalsvm)
#'
#' Fits a survival support vector machine using \pkg{survivalsvm} with configurable
#' loss type, kernel, optimization method, and regularization.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing the variables in the model.
#' @param type Character. One of \code{"regression"}, \code{"vanbelle1"},
#'   \code{"vanbelle2"}, \code{"hybrid"}. Default \code{"regression"}.
#' @param gamma.mu Numeric. Regularization parameter. Default \code{0.1}.
#' @param opt.meth Character. Optimization method; e.g., \code{"quadprog"}.
#' @param kernel Character. Kernel name in \pkg{survivalsvm} (e.g., \code{"lin_kernel"},
#'   \code{"rbf_kernel"}, \code{"add_kernel"}). Default \code{"add_kernel"}.
#' @param diff.meth Optional differentiability method; see \pkg{survivalsvm} docs.
#'
#' @return An \code{"mlsurv_model"} list with elements:
#' \itemize{
#'   \item \code{model} - fitted \code{survivalsvm} object
#'   \item \code{learner} - \code{"survsvm"}
#'   \item \code{engine} - \code{"survivalsvm"}
#'   \item \code{formula}, \code{data}, \code{time}, \code{status}
#' }
#'
#' @details
#' This wrapper standardizes \pkg{survivalsvm} models to the \code{mlsurv_model}
#' interface used across the package.
#'
#' @seealso [predict_survsvm()], \code{survivalsvm::survivalsvm()}
#'
#' @examples
#' \donttest{
#'   fit <- fit_survsvm(
#'     Surv(time, status) ~ age + celltype + karno,
#'     data = veteran, type = "regression", gamma.mu = 0.1, kernel = "lin_kernel"
#'   )
#'   head(predict_survsvm(fit, newdata = veteran[1:5, ], times = c(100, 300, 500)))
#' }
#' @export

fit_survsvm <- function(formula, data,
  type = "regression",
  gamma.mu = 0.1,
  opt.meth = "quadprog",
  kernel = "add_kernel",
  diff.meth = NULL) {

stopifnot(requireNamespace("survivalsvm", quietly = TRUE))
model <- survivalsvm::survivalsvm(
formula = formula,
data = data,
type = type,
gamma.mu = gamma.mu,
opt.meth = opt.meth,
kernel = kernel,
diff.meth = diff.meth)

structure(list(
model = model,
learner = "survsvm",
formula = formula,
data = data,
time = all.vars(formula)[[2]],
status = all.vars(formula)[[3]]
), class = "mlsurv_model", engine = "survivalsvm")
}


#' Predict Survival from a survivalsvm Model
#'
#' Produces survival probabilities at specified time points by converting
#' predicted survival times from \pkg{survivalsvm} using a parametric survival
#' assumption (exponential or Weibull).
#'
#' @param object An \code{mlsurv_model} created by [fit_survsvm()].
#' @param newdata Data frame of new observations.
#' @param times Numeric vector of time points for prediction.
#' @param dist Character. Parametric distribution for mapping predicted times to
#'   survival probabilities: \code{"exp"} (exponential) or \code{"weibull"}.
#'   Default \code{"exp"}.
#' @param shape Numeric shape parameter for Weibull mapping (ignored if \code{dist = "exp"}).
#'   Default \code{1}.
#'
#' @return A \code{data.frame} of survival probabilities with one row per observation
#'   and columns named \code{"t=<time>"}.
#'
#' @details
#' \pkg{survivalsvm} returns predicted survival times. This wrapper converts them
#' to survival probabilities at \code{times} via:
#' \itemize{
#'   \item Exponential: \eqn{S(t) = \exp(-t/\hat{T})}
#'   \item Weibull: \eqn{S(t) = \exp\{-(t/\hat{T})^{\text{shape}}\}}
#' }
#' where \eqn{\hat{T}} is the predicted time. Choose \code{dist} and \code{shape}
#' to match your application.
#'
#' @seealso [fit_survsvm()], \code{survivalsvm::predict.survivalsvm()}
#'
#' @examples
#' \donttest{
#'   fit <- fit_survsvm(
#'     Surv(time, status) ~ age + celltype + karno,
#'     data = veteran, kernel = "lin_kernel"
#'   )
#'   times <- c(100, 300, 500)
#'   predict_survsvm(fit, newdata = veteran[1:5, ], times = times, dist = "exp")
#'   predict_survsvm(fit, newdata = veteran[1:5, ], times = times, dist = "weibull", shape = 1.5)
#' }
#' @export

predict_survsvm <- function(object, newdata, times, dist = "exp", shape = 1) {

if (!is.null(object$learner) && object$learner != "survsvm") {
warning("Object passed to predict_survsvm() may not come from fit_survsvm().")
}

stopifnot(requireNamespace("survivalsvm", quietly = TRUE))

newdata <- as.data.frame(newdata)
pred <- predict(object$model, newdata = newdata)
predicted_times <- as.numeric(pred$predicted)

dist <- match.arg(dist, choices = c("exp", "weibull"))

survmat <- switch(dist,
"exp" = sapply(times, function(t) exp(-t / predicted_times)),
"weibull" = sapply(times, function(t) exp(- (t / predicted_times)^shape))
)

colnames(survmat) <- paste0("t=", times)
as.data.frame(survmat)
}


#' Hyperparameter Tuning for survivalsvm
#'
#' Cross-validated tuning over a grid of \pkg{survivalsvm} hyperparameters.
#' Returns a results table or refits and returns the best model.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data Training data frame.
#' @param times Numeric vector of time points used for evaluation.
#' @param metrics Character vector of metric names to aggregate (e.g., \code{c("cindex","ibs")}).
#' @param param_grid Named list of candidate values. Example:
#'   \code{list(gamma.mu = c(0.01, 0.1), kernel = c("lin_kernel","add_kernel"))}.
#' @param folds Integer. Number of CV folds. Default \code{5}.
#' @param seed Integer. Random seed for reproducibility.
#' @param refit_best Logical. If \code{TRUE}, refit best config on full data and return it.
#' @param dist Character. Mapping distribution passed to [predict_survsvm()] (\code{"exp"} or \code{"weibull"}).
#' @param shape Numeric. Weibull shape for mapping (used if \code{dist = "weibull"}).
#'
#' @return
#' If \code{refit_best = FALSE}, a \code{tibble} with one row per configuration and
#' columns for aggregated \code{metrics}. If \code{refit_best = TRUE}, an
#' \code{mlsurv_model} with attribute \code{"tuning_results"} containing the full table.
#'
#' @details
#' For stability, the tuner fixes \code{type = "regression"}, \code{opt.meth = "quadprog"},
#' and \code{diff.meth = NULL}. The prediction mapping is controlled by \code{dist} and \code{shape}.
#'
#' @seealso [fit_survsvm()], [predict_survsvm()], [cv_survlearner()]
#'
#' @examples
#' \donttest{
#'   grid <- list(
#'     gamma.mu = c(0.01, 0.1),
#'     kernel = c("lin_kernel", "add_kernel")
#'   )
#'   # Return best refitted model
#'   best <- tune_survsvm(
#'     formula = Surv(time, status) ~ age + celltype + karno,
#'     data = veteran,
#'     times = c(100, 300, 500),
#'     metrics = c("cindex", "ibs"),
#'     param_grid = grid,
#'     folds = 3,
#'     refit_best = TRUE
#'   )
#'   predict_survsvm(best, newdata = veteran[1:5, ], times = c(100,300))
#'
#'   # Or get the full results table
#'   res <- tune_survsvm(
#'     formula = Surv(time, status) ~ age + celltype + karno,
#'     data = veteran,
#'     times = c(100, 300, 500),
#'     metrics = c("cindex", "ibs"),
#'     param_grid = grid,
#'     folds = 3,
#'     refit_best = FALSE
#'   )
#'   res
#' }
#' @export

tune_survsvm <- function(formula, data, times,
   metrics = "cindex",
   param_grid,
   folds = 5, seed = 42,
   refit_best = FALSE,
   dist = "exp", shape = 1) {
stopifnot(is.list(param_grid))
param_df <- tidyr::crossing(!!!param_grid)

# Fixed safe defaults
fixed_type <- "regression"
fixed_opt_meth <- "quadprog"
fixed_diff_meth <- NULL

results <- purrr::pmap_dfr(param_df, function(gamma.mu, kernel) {
set.seed(seed)
tryCatch({
res_cv <- cv_survlearner(
formula = formula,
data = data,
fit_fun = fit_survsvm,
pred_fun = function(...) predict_survsvm(..., dist = dist, shape = shape),
times = times,
metrics = metrics,
folds = folds,
seed = seed,
gamma.mu = gamma.mu,
kernel = kernel,
type = fixed_type,
opt.meth = fixed_opt_meth,
diff.meth = fixed_diff_meth
)

summary <- cv_summary(res_cv)

tibble::tibble(
gamma.mu = gamma.mu,
kernel = kernel
) |>
dplyr::bind_cols(
tidyr::pivot_wider(summary[, c("metric", "mean")],
   names_from = metric, values_from = mean)
)
}, error = function(e) {
invisible(cat(glue::glue("Skipping (gamma.mu={gamma.mu}, kernel={kernel}): {e$message}\n"),
file = "errors.log", append = TRUE)) ## keep a trace without showing it on the console
NULL
})
})

if (nrow(results) == 0) {
warning("All tuning combinations failed.")
return(NULL)
}

results <- dplyr::arrange(results, dplyr::desc(.data[[metrics[1]]]))

if (refit_best) {
best <- results[1, ]
model <- fit_survsvm(
formula = formula,
data = data,
gamma.mu = best$gamma.mu,
kernel = best$kernel,
type = fixed_type,
opt.meth = fixed_opt_meth,
diff.meth = fixed_diff_meth
)
attr(model, "tuning_results") <- results
class(model) <- c("tune_surv", class(model))
return(model)
}

class(results) <- c("tune_surv", class(results))
return(results)
}


