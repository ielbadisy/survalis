
#' Concordance Index from a Survival-Probability Matrix
#'
#' Computes Harrell's concordance index using predicted survival probabilities
#' at a specified time point (or the last column if \code{t_star} is \code{NULL}).
#'
#' @param object A \code{\link[survival]{Surv}} object of length \eqn{n}.
#' @param predicted An \code{n x k} matrix or data frame of survival probabilities
#'   with columns named \code{"t=<time>"}.
#' @param t_star Optional numeric time at which to evaluate the c-index; if omitted,
#'   the rightmost column of \code{predicted} is used.
#'
#' @details
#' Risk scores are defined as \code{1 - S(t)} at the chosen time. Ties receive
#' partial credit (0.5). Pairs not comparable due to censoring are excluded.
#'
#' @return A named numeric scalar: \code{"C index"}.
#'
#' @examples
#' # cindex_survmat(Surv(time, status), sp_matrix, t_star = 365)
#' @export

cindex_survmat <- function(object, predicted, t_star = NULL) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")
  time <- object[, 1]
  status <- object[, 2]

  # Extract column for the given time point
  if (!is.null(t_star)) {
    t_name <- paste0("t=", t_star)
    if (!(t_name %in% colnames(predicted))) {
      stop("t_star = ", t_star, " not found in predicted survival matrix.")
    }
    surv_prob <- if (is.matrix(predicted)) predicted[, t_name] else predicted[[t_name]]
  } else {
    surv_prob <- if (is.matrix(predicted)) predicted[, ncol(predicted)] else predicted[[ncol(predicted)]]
  }

  risk_score <- 1 - surv_prob
  permissible <- 0
  concord <- 0
  par_concord <- 0
  n <- length(time)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((time[i] < time[j] & status[i] == 0) | (time[j] < time[i] & status[j] == 0)) next
      if (time[i] == time[j] & status[i] == 0 & status[j] == 0) next
      permissible <- permissible + 1
      if (time[i] != time[j]) {
        if ((time[i] < time[j] & risk_score[i] > risk_score[j]) |
            (time[j] < time[i] & risk_score[j] > risk_score[i])) {
          concord <- concord + 1
        } else if (risk_score[i] == risk_score[j]) {
          par_concord <- par_concord + 0.5
        }
      } else {
        if (status[i] + status[j] > 0) {
          if (risk_score[i] == risk_score[j]) concord <- concord + 1
          else par_concord <- par_concord + 0.5
        }
      }
    }
  }

  C_index <- (concord + par_concord) / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}

#' Brier Score with IPCW for a Single Time Point
#'
#' Computes the inverse-probability-of-censoring weighted (IPCW) Brier score
#' at a single time \code{t_star}.
#'
#' @param object A \code{\link[survival]{Surv}} object.
#' @param pre_sp Numeric vector of predicted survival probabilities \eqn{S(t^{*}\mid x_i)}.
#' @param t_star Numeric evaluation time.
#'
#' @details
#' The censoring distribution \eqn{G(t)} is estimated via Kaplan-Meier on
#' \code{1 - status}. Observed events before \code{t_star} contribute
#' \eqn{S(t_i)^2 / G(t_i)}; those at risk at \code{t_star} contribute
#' \eqn{(1 - S(t^{*}))^2 / G(t^{*})}. Returns \code{NA} if \eqn{G(t^{*})} is
#' undefined or zero.
#'
#' @return A named numeric scalar: \code{"brier"}.
#'
#' @examples
#' # brier(Surv(time, status), pre_sp = sp[,1], t_star = 180)
#' @export

brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  time <- object[, 1]
  status <- object[, 2]
  km_fit <- survfit(Surv(time, 1 - status) ~ 1)

  all_times <- sort(unique(c(t_star, time)))
  G_all <- summary(km_fit, times = all_times, extend = TRUE)$surv
  names(G_all) <- as.character(all_times)

  Gt <- G_all[as.character(t_star)]

  if (is.na(Gt) || Gt == 0) {
    warning("Brier: Gt(t_star) is NA or zero at t = ", t_star, "; returning NA.")
    return(NA_real_)
  }

  early_event <- which(time < t_star & status == 1)
  later_risk  <- which(time >= t_star)
  score_vec <- numeric(length(time))

  if (length(early_event) > 0) {
    Gti <- G_all[as.character(time[early_event])]
    valid <- which(!is.na(Gti) & Gti > 0)
    score_vec[early_event[valid]] <- (pre_sp[early_event[valid]]^2) / Gti[valid]
  }

  if (Gt > 0 && length(later_risk) > 0) {
    score_vec[later_risk] <- ((1 - pre_sp[later_risk])^2) / Gt
  }

  bs_value <- mean(score_vec, na.rm = TRUE)
  names(bs_value) <- "brier"
  return(round(bs_value, 6))
}

#' Integrated Brier Score (Discrete Integration)
#'
#' Computes the Integrated Brier Score (IBS) over a vector of times, using
#' discrete left Riemann integration of IPCW Brier scores.
#'
#' @param object A \code{\link[survival]{Surv}} object.
#' @param sp_matrix Matrix/data frame of survival probabilities with one column
#'   per time in \code{times}.
#' @param times Numeric vector of strictly increasing times (length must equal
#'   \code{ncol(sp_matrix)}).
#'
#' @return A named numeric scalar: \code{"ibs"}.
#'
#' @examples
#' # ibs_survmat(Surv(time, status), sp_matrix = sp, times = c(90,180,365))
#' @export

ibs_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")
  if (length(times) == 1) {
    return(brier(object, pre_sp = sp_matrix[, 1], t_star = times[1]))
  }

  times <- sort(times)
  brier_vec <- vapply(seq_along(times), function(i) {
    brier(object, pre_sp = sp_matrix[, i], t_star = times[i])
  }, numeric(1))

  dt <- diff(times)
  integral <- sum(brier_vec[-length(brier_vec)] * dt)
  ibs_value <- integral / (max(times) - min(times))
  names(ibs_value) <- "ibs"
  return(round(ibs_value, 6))
}

#' Integrated Absolute Error Against Kaplan-Meier
#'
#' Computes \eqn{\int \mid \bar{S}(t) - \hat{S}_{KM}(t) \mid \, dt} where
#' \eqn{\bar{S}(t)} is the mean predicted survival across subjects and
#' \eqn{\hat{S}_{KM}(t)} is the Kaplan-Meier estimate. Integration is carried
#' out over KM event times using a left Riemann sum.
#'
#' @param object A \code{\link[survival]{Surv}} object.
#' @param sp_matrix Matrix/data frame of survival probabilities (rows = subjects,
#'   columns aligned with \code{times}).
#' @param times Numeric vector of times corresponding to the columns of \code{sp_matrix}.
#'
#' @return A named numeric scalar: \code{"iae"}.
#'
#' @examples
#' # iae_survmat(Surv(time, status), sp, times)
#' @export

iae_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  pred_at_km <- approx(x = times, y = mean_pred_surv, xout = km_data$time,
                       method = "constant", rule = 2)$y

  dt <- diff(km_data$time)
  t_IAE <- abs(km_data$surv - pred_at_km)
  iae <- sum(t_IAE[-length(t_IAE)] * dt)
  names(iae) <- "iae"
  return(round(iae, 4))
}

#' Integrated Squared Error Against Kaplan-Meier
#'
#' Computes \eqn{\int ( \bar{S}(t) - \hat{S}_{KM}(t) )^2 \, dt} where
#' \eqn{\bar{S}(t)} is the mean predicted survival and \eqn{\hat{S}_{KM}(t)}
#' is the Kaplan-Meier curve. Integration uses KM event times and a left
#' Riemann sum.
#'
#' @param object A \code{\link[survival]{Surv}} object.
#' @param sp_matrix Matrix/data frame of survival probabilities (rows = subjects,
#'   columns aligned with \code{times}).
#' @param times Numeric vector of times corresponding to the columns of \code{sp_matrix}.
#'
#' @return A named numeric scalar: \code{"ise"}.
#'
#' @examples
#' # ise_survmat(Surv(time, status), sp, times)
#' @export

ise_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  pred_at_km <- approx(x = times, y = mean_pred_surv, xout = km_data$time,
                       method = "constant", rule = 2)$y

  dt <- diff(km_data$time)
  t_ISE <- (km_data$surv - pred_at_km)^2
  ise <- sum(t_ISE[-length(t_ISE)] * dt)
  names(ise) <- "ise"
  return(round(ise, 4))
}


#' Cross-Validate a Survival Learner (fold-mapped with \code{fmapn})
#'
#' Runs k-fold cross-validation for any pair of \code{fit_fun}/\code{pred_fun}
#' that follow the package's learner contracts, and returns tidy per-fold metric values.
#' Fold iteration is handled by \code{functionals::fmapn()} with optional parallel
#' execution and a progress bar.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing all variables in \code{formula}.
#' @param fit_fun Function with signature \code{fit_fun(formula, data, ...)} that
#'   returns an \code{mlsurv_model}.
#' @param pred_fun Function with signature \code{pred_fun(object, newdata, times, ...)}
#'   returning a survival-probability matrix/data frame (columns named \code{"t=<time>"}).
#' @param times Numeric vector of evaluation times (passed to \code{pred_fun}).
#' @param metrics Character vector of metrics to compute. Supported:
#'   \code{"cindex"}, \code{"brier"} (single time), \code{"ibs"}, \code{"iae"}, \code{"ise"}.
#' @param folds Integer; number of folds (default \code{5}).
#' @param seed Integer random seed for reproducibility (default \code{123}).
#' @param verbose Logical; print row-dropping due to missingness (default \code{FALSE}).
#' @param ncores Integer; number of CPU cores for \code{fmapn} parallel mapping
#'   (default \code{1}). Set \code{> 1} to enable parallelism.
#' @param pb Logical; show a progress bar during fold mapping (default \code{TRUE}).
#' @param ... Additional arguments forwarded to \code{fit_fun}.
#'
#' @details
#' The routine:
#' \enumerate{
#' \item Validates \code{Surv(...)} on the LHS and warns against using \code{.} in formulas.
#' \item Drops rows with missing values in any variables referenced by \code{formula}.
#' \item Supports \code{Surv(time, status == k)} by recoding the status to 0/1.
#' \item Builds stratified v-folds on the status indicator (\pkg{rsample}).
#' \item For each fold: fits on the analysis set, predicts on the assessment set, and computes metrics.
#' }
#'
#' Fold iteration is performed via \code{functionals::fmapn()}, which preserves
#' per-fold identifiers (\code{id}, \code{fold}) and returns a list ready for
#' \code{dplyr::bind_rows()}.
#'
#' @return A tibble with columns: \code{splits} (rsample split object),
#'   \code{id}, \code{fold}, \code{metric}, and \code{value}.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(veteran)
#' cv_results <- cv_survlearner(
#'   formula = Surv(time, status) ~ karno,
#'   data = veteran,
#'   fit_fun = fit_bart,
#'   pred_fun = predict_bart,
#'   times = c(10, 100, 600, 900),
#'   metrics = c("cindex", "ibs"),
#'   folds = 5,
#'   ncores = max(1, parallel::detectCores() - 1),
#'   pb = TRUE
#' )
#' }
#'
#' @seealso \code{\link[functionals]{fmapn}}
#' @importFrom functionals fmapn
#' @keywords survival cross-validation fmapn parallel
#' @export



cv_survlearner <- function(formula, data,
  fit_fun, pred_fun,
  times,
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 123,
  verbose = FALSE,
  ncores = 1,
  pb = TRUE,
  ...) {

if ("." %in% all.vars(update(formula, . ~ 0))) {
warning("Please avoid using '.' in the formula. Specify all predictors explicitly.")
}

tf <- terms(formula, data = data)
outcome <- attr(tf, "variables")[[2]]

if (!inherits(outcome, "call") || outcome[[1]] != as.name("Surv")) {
stop("Formula must begin with a Surv(...) outcome.")
}

time_col <- as.character(outcome[[2]])
status_expr <- outcome[[3]]

if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
status_col <- as.character(status_expr[[2]])
event_value <- eval(status_expr[[3]], data)
recode_status <- TRUE
} else {
status_col <- as.character(status_expr)
event_value <- 1
recode_status <- FALSE
}

all_vars <- all.vars(formula)
n_before <- nrow(data)
data <- tidyr::drop_na(data, dplyr::all_of(all_vars))
n_after <- nrow(data)

if (verbose && n_after < n_before) {
message("Dropped ", n_before - n_after, " rows with missing values (from ", n_before, " to ", n_after, ").")
}

if (recode_status) {
data[[status_col]] <- as.integer(data[[status_col]] == event_value)
}

set.seed(seed)
folds_data <- rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col))

# cross-validation loop via fmapn
results <- functionals::fmapn(
list(split = folds_data$splits, id = folds_data$id, fold = seq_len(folds)),
function(split, id, fold) {
train <- rsample::analysis(split)
test  <- rsample::assessment(split)

model <- fit_fun(formula = formula, data = train, ...)
pred  <- pred_fun(model, newdata = test, times = times)

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
dplyr::mutate(value = lapply(metric, function(metric) {
switch(metric,
"cindex" = cindex_survmat(surv_obj, predicted = pred, t_star = max(times)),
"brier"  = {
if (length(times) != 1) stop("Brier requires a single time point.")
brier(surv_obj, pre_sp = pred[, 1], t_star = times)
},
"ibs" = ibs_survmat(surv_obj, sp_matrix = pred, times = times),
"iae" = iae_survmat(surv_obj, sp_matrix = pred, times = times),
"ise" = ise_survmat(surv_obj, sp_matrix = pred, times = times),
stop("Unknown metric: ", metric)
)
})) |>
tidyr::unnest(cols = value) |>
dplyr::mutate(splits = list(split), id = id, fold = fold) |>
dplyr::relocate(splits, id, fold)
},
ncores = ncores, pb = pb
)
  dplyr::bind_rows(results)
}



#' Summarize Cross-Validation Results
#'
#' Produces mean, standard deviation, standard error, and 95\% Wald intervals
#' for each metric returned by \code{\link{cv_survlearner}}.
#'
#' @param cv_results A tibble/data frame as returned by \code{cv_survlearner()}.
#'
#' @return A tibble with columns: \code{metric}, \code{mean}, \code{sd}, \code{n},
#'   \code{se}, \code{lower}, \code{upper}.
#'
#' @examples
#' # cv_summary(cv_results)
#' @export

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


#' Boxplot of Cross-Validation Metric Distributions
#'
#' Visualizes per-fold metric values from \code{\link{cv_survlearner}} using a
#' boxplot with jittered points.
#'
#' @param cv_results A tibble/data frame as returned by \code{cv_survlearner()}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' # cv_plot(cv_results)
#' @export

cv_plot <- function(cv_results) {
  ggplot(cv_results, aes(x = metric, y = value)) +
    geom_boxplot(fill = "skyblue", outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    labs(title = "Cross-Validation Performance", x = "Metric", y = "Value") +
    theme_minimal()
}


#' Score a Fitted Survival Model on Its Training Data
#'
#' Computes one or more performance metrics for a fitted \code{mlsurv_model},
#' predicting on the training data with the model's corresponding \code{predict_*}
#' function.
#'
#' @param model An object of class \code{"mlsurv_model"}.
#' @param times Numeric vector of evaluation times. For \code{"brier"}, must be a single time.
#' @param metrics Character vector of metrics to compute. Supported:
#'   \code{"cindex"}, \code{"ibs"}, \code{"brier"}, \code{"iae"}, \code{"ise"}.
#'
#' @details
#' The function constructs the appropriate \code{predict_*} function name from
#' \code{model$learner}, predicts survival probabilities on \code{model$data},
#' builds a \code{Surv} object from \code{model$formula}, and computes the metrics.
#' If \code{"brier"} is requested with multiple \code{times}, an error is thrown.
#'
#' @return A tibble with columns \code{metric} and \code{value}.
#'
#' @examples
#' # score_survmodel(fitted_model, times = c(90, 180), metrics = c("ibs", "cindex"))
#' @export

score_survmodel <- function(model, times, metrics = c("cindex", "ibs", "brier", "iae", "ise")) {
  stopifnot(inherits(model, "mlsurv_model"))
  stopifnot(!missing(times))

  # check brier with multiple times
  if ("brier" %in% metrics && length(times) != 1) {
    stop("The Brier score requires a single evaluation time.\n",
         "Please use `times = t` with a single value or remove 'brier' from the metric list.")
  }

  # extract model fields
  learner <- model$learner
  formula <- model$formula
  data    <- model$data

  if (is.null(learner)) stop("Model object must have a 'learner' field.")

  # get prediction function
  pred_fun_name <- paste0("predict_", learner)
  if (!exists(pred_fun_name, mode = "function")) {
    stop("Prediction function not found: ", pred_fun_name)
  }
  pred_fun <- get(pred_fun_name, mode = "function")
  sp_matrix <- pred_fun(model, newdata = data, times = times)

  # extract Surv object
  tf <- terms(formula, data = data)
  outcome <- attr(tf, "variables")[[2]]
  time_col <- as.character(outcome[[2]])
  status_expr <- outcome[[3]]

  if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
    status_col <- as.character(status_expr[[2]])
    event_value <- eval(status_expr[[3]], data)
    status_vector <- as.integer(data[[status_col]] == event_value)
  } else {
    status_col <- as.character(status_expr)
    event_value <- 1
    status_vector <- data[[status_col]]
  }

  surv_obj <- survival::Surv(time = data[[time_col]], event = status_vector)

  # score each metric
  tibble::tibble(metric = metrics) |>
    dplyr::mutate(value = purrr::map(metric, function(metric) {
      switch(metric,
             "cindex" = cindex_survmat(surv_obj, predicted = sp_matrix, t_star = max(times)),
             "brier"  = brier(surv_obj, pre_sp = sp_matrix[, 1], t_star = times),
             "ibs"    = ibs_survmat(surv_obj, sp_matrix, times),
             "iae"    = iae_survmat(surv_obj, sp_matrix, times),
             "ise"    = ise_survmat(surv_obj, sp_matrix, times),
             stop("Unknown metric: ", metric)
      )
    })) |>

    tidyr::unnest(cols = value)
  }



