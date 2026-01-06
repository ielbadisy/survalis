#' Fit a Stacked Survival Meta-Learner (Time-Varying NNLS)
#'
#' Learns nonnegative, time-specific stacking weights across a set of base
#' survival learners using a nonnegative least squares (NNLS) approach.
#' For each evaluation time \eqn{t^{*}}, the target is the survival indicator
#' \eqn{Y = I(T > t^{*})}, and the predictors are the base learners’ survival
#' probabilities \eqn{S_\ell(t^{*} \mid x)}. Weights are constrained to be
#' nonnegative and renormalized to sum to 1 per time point.
#'
#' @param base_preds A named list of matrices or data frames, one per base learner,
#'   each of dimension \code{n × length(times)}, with columns named
#'   \code{"t=<time>"} corresponding to evaluation times.
#' @param time Numeric vector of observed survival times (length \code{n}).
#' @param status Binary event indicator (1 = event, 0 = censored) of length \code{n}.
#' @param times Numeric vector of time points at which to learn stacking weights.
#' @param base_models A named list of fitted base learner models. Names must match
#'   those in \code{base_preds} and correspond to available \code{predict_<learner>()}
#'   functions.
#' @param formula A \code{Surv()} formula used to fit the base models; stored for metadata
#'   and downstream scoring.
#' @param data The training data frame used to fit the base learners.
#'
#' @details
#' For each time point \eqn{t \in \texttt{times}}, this function solves
#' \deqn{
#'   \min_{w \ge 0} \|Y - S w\|_2^2
#'   \quad \text{s.t. } \sum_j w_j = 1,
#' }
#' where \eqn{Y = I(T > t)} and \eqn{S} is the matrix of base learner survival
#' probabilities at time \eqn{t}. The nonnegative least squares problem is solved
#' using \pkg{nnls}. If the NNLS solution sums to zero or contains non-finite values,
#' uniform weights (\eqn{1/L}) are assigned across learners.
#'
#' This approach produces a *time-varying ensemble* that adaptively reweights the
#' contributions of each base learner across the survival horizon.
#'
#' @return An object of class \code{c("mlsurv_model","survmetalearner")} containing:
#' \describe{
#'   \item{weights}{Matrix (\eqn{L × T}) of nonnegative weights,
#'     rows = base learners, columns = \code{"t=<time>"} points.}
#'   \item{base_models}{Named list of base learner fits.}
#'   \item{base_preds}{Named list of base prediction matrices on the training set.}
#'   \item{learners}{Character vector of base learner names.}
#'   \item{formula, data, time, status}{Training metadata for reproducibility.}
#'   \item{learner}{String \code{"survmetalearner"} used for dispatch in
#'     \code{predict_survmetalearner()}.}
#' }
#'
#' @section Implementation Notes:
#' * Columns in \code{base_preds} must exactly match \code{paste0("t=", times)}.
#' * If fewer than two base learners are provided, the function falls back to
#'   uniform weighting.
#' * All base learners are assumed to produce calibrated survival probabilities.
#'
#' @seealso
#' \code{\link{predict_survmetalearner}},
#' \code{\link{plot_survmetalearner_weights}},
#' \code{\link{cv_survmetalearner}}
#'
#' @examples
#' \donttest{
#' veteran < survival::veteran
#' veteran$status <- as.integer(veteran$status == 1)
#' form <- Surv(time, status) ~ age + karno + celltype
#' times <- c(100, 200, 400)
#'
#' mod1 <- fit_coxph(form, data = veteran)
#' mod2 <- fit_rpart(form, data = veteran)
#'
#' base_models <- list(coxph = mod1, rpart = mod2)
#' base_preds <- list(
#'   coxph = predict_coxph(mod1, veteran, times),
#'   rpart = predict_rpart(mod2, veteran, times)
#' )
#'
#' meta <- fit_survmetalearner(
#'   base_preds = base_preds, time = veteran$time, status = veteran$status,
#'   times = times, base_models = base_models, formula = form, data = veteran
#' )
#' meta$weights
#' }
#'
#' @export

fit_survmetalearner <- function(base_preds, time, status, times,
                                base_models, formula, data) {
  stopifnot(is.list(base_preds), is.list(base_models), is.data.frame(data))
  if (!setequal(names(base_preds), names(base_models))) {
    stop("names(base_preds) must match names(base_models).")
  }

  # ensure matrices and column alignment to requested times
  tnames <- paste0("t=", times)
  base_preds <- lapply(base_preds, function(x) {
    M <- if (!is.matrix(x)) as.matrix(x) else x
    if (!all(tnames %in% colnames(M))) {
      stop("Each base_preds matrix must have columns named exactly ", paste(tnames, collapse = ", "))
    }
    M[, tnames, drop = FALSE]
  })

  learners <- names(base_preds)
  T <- length(times)
  weight_matrix <- matrix(NA_real_, nrow = length(learners), ncol = T,
                          dimnames = list(learners, tnames))

  for (j in seq_along(times)) {
    t_star <- times[j]
    Y <- as.integer(time > t_star)
    X <- sapply(base_preds, function(S) S[, j])
    fit <- nnls::nnls(as.matrix(X), Y)
    w <- coef(fit)
    s <- sum(w)
    w <- if (!is.finite(s) || s <= 0) rep(1 / length(w), length(w)) else w / s
    weight_matrix[, j] <- w
  }

  structure(list(
    weights = weight_matrix,
    base_models = base_models,
    base_preds = base_preds,
    formula = formula,
    data = data,
    time = time,
    status = status,
    learner = "survmetalearner",
    learners = learners
  ), class = c("mlsurv_model", "survmetalearner"), engine = "survalis")
}


#' Predict with a Stacked Survival Meta-Learner
#'
#' Generates stacked survival probabilities by combining base learner predictions
#' using the time-varying nonnegative weights learned via
#' \code{\link{fit_survmetalearner}}.
#'
#' @param model A fitted \code{"survmetalearner"} object produced by
#'   \code{\link{fit_survmetalearner}}.
#' @param newdata A data frame of new observations for prediction.
#' @param times Numeric vector of evaluation times at which to compute survival
#'   probabilities. Must be a subset of the training times used to fit the
#'   meta-learner.
#'
#' @details
#' For each base learner listed in \code{model$learners}, the corresponding
#' \code{predict_<learner>} function is called to obtain survival probabilities
#' \eqn{S_\ell(t \mid x)} for each requested time point.
#'
#' The meta-prediction at each time \eqn{t} is computed as:
#' \deqn{
#'   \hat{S}(t \mid x) = \sum_{\ell=1}^{L} w_\ell(t) S_\ell(t \mid x),
#' }
#' where \eqn{w_\ell(t)} are the nonnegative weights estimated at training time
#' (see \code{\link{fit_survmetalearner}}).
#'
#' - If any requested time points are not present in the learned weight matrix,
#'   an error is raised.
#' - When only one base learner is available, or if all learned weights at a
#'   given time point are zero or non-finite, uniform weights (\eqn{1/L}) are used.
#' - Each base learner’s prediction function must accept arguments
#'   \code{(object, newdata, times)} and return a numeric matrix/data frame
#'   with columns named \code{"t=<time>"}.
#'
#' Internally, this function safeguards against matrix shape issues and
#' ensures that all predicted probabilities remain bounded in \eqn{[0, 1]}.
#'
#' @return
#' A data frame of survival probabilities with:
#' \itemize{
#'   \item one row per observation in \code{newdata};
#'   \item one column per requested time (columns named \code{"t=<time>"}).
#' }
#'
#' @section Implementation Notes:
#' * Predictions are purely aggregative: the meta-learner performs no re-fitting.
#' * If only one base learner was used, the predictions will exactly match that
#'   learner's output.
#' * This function is compatible with all base learners implementing
#'   \code{predict_<learner>()} returning survival probabilities.
#'
#' @seealso
#' \code{\link{fit_survmetalearner}},
#' \code{\link{plot_survmetalearner_weights}},
#' \code{\link{cv_survmetalearner}}
#'
#' @examples
#' \donttest{
#' veteran < survival::veteran
#' veteran$status <- as.integer(veteran$status == 1)
#' form <- Surv(time, status) ~ age + karno + celltype
#' times <- c(100, 200, 400)
#'
#' # fit base learners
#' mod1 <- fit_coxph(form, data = veteran)
#' mod2 <- fit_rpart(form, data = veteran)
#'
#' base_models <- list(coxph = mod1, rpart = mod2)
#' base_preds  <- list(
#'   coxph = predict_coxph(mod1, veteran, times),
#'   rpart = predict_rpart(mod2, veteran, times)
#' )
#'
#' meta <- fit_survmetalearner(
#'   base_preds = base_preds, time = veteran$time, status = veteran$status,
#'   times = times, base_models = base_models, formula = form, data = veteran
#' )
#'
#' # generate predictions on new data
#' preds <- predict_survmetalearner(meta, newdata = veteran[1:5, ], times = c(100, 200))
#' head(preds)
#' }
#'
#' @export

predict_survmetalearner <- function(model, newdata, times) {
  stopifnot(inherits(model, "survmetalearner"))

  learners <- model$learners
  base_models <- model$base_models
  requested_tnames <- paste0("t=", times)

  # requested times must have weights!!
  if (!all(requested_tnames %in% colnames(model$weights))) {
    stop("Requested 'times' must be a subset of the training times used by the meta-learner.")
  }
  W <- model$weights[, requested_tnames, drop = FALSE]

  # get base predictions for each learner
  base_preds_test <- lapply(learners, function(learner) {
    pred_fun <- get(paste0("predict_", learner), mode = "function")
    pred <- pred_fun(base_models[[learner]], newdata = newdata, times = times)
    M <- as.matrix(pred)
    colnames(M) <- requested_tnames
    M
  })
  names(base_preds_test) <- learners

  n <- nrow(base_preds_test[[1]])
  T <- length(times)
  survmat <- matrix(NA_real_, nrow = n, ncol = T)

  for (j in seq_len(T)) {
    # robust n*L base matrix at time j (no drop when L=1, this fixed the single base learner case issue)  
    base_mat_j <- do.call(cbind, lapply(base_preds_test, function(M) M[, j, drop = TRUE]))
    if (is.null(dim(base_mat_j))) base_mat_j <- cbind(base_mat_j)

    wj <- W[, j]
    if (!all(is.finite(wj)) || sum(wj) <= 0) wj <- rep(1 / length(wj), length(wj))

    survmat[, j] <- as.vector(base_mat_j %*% wj)
  }

  colnames(survmat) <- requested_tnames
  as.data.frame(survmat)
}


#' Plot Time-Varying Stacking Weights from a Survival Meta-Learner
#'
#' Visualizes how each base learner’s contribution evolves over time, as learned
#' by the time-varying NNLS stacking model from
#' \code{\link{fit_survmetalearner}}.
#'
#' @param model A fitted \code{"survmetalearner"} object produced by
#'   \code{\link{fit_survmetalearner}} or returned within
#'   \code{\link{cv_survmetalearner}}.
#'
#' @details
#' This plot displays the nonnegative ensemble weights \eqn{w_\ell(t)} assigned
#' to each base learner \eqn{\ell} across the evaluation times \eqn{t} used in
#' training.
#'
#' Each curve corresponds to one learner, showing how its relative importance
#' changes over the survival horizon. The weights at each \eqn{t} sum to 1,
#' reflecting the convex NNLS constraint in the meta-learner.
#'
#' Internally, the function reshapes the weight matrix (\eqn{L × T}) to long
#' format and draws smooth trajectories with \pkg{ggplot2}.  
#'
#' This visualization is useful to:
#' \itemize{
#'   \item interpret dynamic model weighting behavior;
#'   \item identify learners that dominate at early vs. late follow-up;
#'   \item assess stability and complementarity among base learners.
#' }
#'
#' @return
#' A \pkg{ggplot2} object displaying the time-varying weights, with:
#' \itemize{
#'   \item \emph{x-axis}: time points (\code{as.numeric(times)});
#'   \item \emph{y-axis}: normalized NNLS weights;
#'   \item \emph{color}: base learner identity.
#' }
#'
#' @seealso
#' \code{\link{fit_survmetalearner}},
#' \code{\link{predict_survmetalearner}},
#' \code{\link{cv_survmetalearner}}
#'
#' @examples
#' \donttest{
#' veteran < survival::veteran
#' veteran$status <- as.integer(veteran$status == 1)
#' form <- Surv(time, status) ~ age + karno + celltype
#' times <- c(100, 200, 400)
#'
#' mod1 <- fit_coxph(form, data = veteran)
#' mod2 <- fit_rpart(form, data = veteran)
#'
#' base_models <- list(coxph = mod1, rpart = mod2)
#' base_preds  <- list(
#'   coxph = predict_coxph(mod1, veteran, times),
#'   rpart = predict_rpart(mod2, veteran, times)
#' )
#'
#' meta <- fit_survmetalearner(
#'   base_preds = base_preds, time = veteran$time, status = veteran$status,
#'   times = times, base_models = base_models, formula = form, data = veteran
#' )
#'
#' plot_survmetalearner_weights(meta)
#' }
#'
#' @export

plot_survmetalearner_weights <- function(model) {
  stopifnot(inherits(model, "survmetalearner"))

  W <- as.data.frame(model$weights)
  W$learner <- rownames(model$weights)

  W_long <- tidyr::pivot_longer(W, -learner, names_to = "time", values_to = "weight") |>
    dplyr::mutate(time = as.numeric(sub("t=", "", time)))

  ggplot2::ggplot(W_long, ggplot2::aes(x = time, y = weight, color = learner)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::labs(x = "Time", y = "Weight", title = "NNLS stacking weights over time") +
    ggplot2::theme_minimal(base_size = 14)
}


#' Cross-Validate a Stacked Survival Meta-Learner
#'
#' Performs v-fold cross-validation for the NNLS stacking meta-learner over a fixed
#' set of base learners and evaluation times. On each fold, time-varying weights
#' are learned on the analysis set and evaluated on the assessment set.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data Data frame containing variables in \code{formula}.
#' @param times Numeric vector of evaluation times for stacking and scoring.
#' @param base_models Named list of fitted base learner models. Names must match
#'   \code{names(base_preds)} and correspond to available \code{predict_<learner>()}.
#' @param base_preds Named list of training-set prediction matrices/data frames with
#'   columns \code{"t=<time>"} matching \code{times} and rows aligned with \code{data}.
#' @param folds Integer; number of v-folds (default \code{5}).
#' @param metrics Character vector of metrics to compute. Supported:
#'   \code{"cindex"}, \code{"brier"} (single time), \code{"ibs"}, \code{"iae"}, \code{"ise"}.
#' @param seed Integer random seed (default \code{123}).
#' @param verbose Logical; print fold progress (default \code{TRUE}).
#' @param ncores Integer; number of cores for fold-level parallelism via
#'   \code{functionals::fmapn()} (default \code{1} = sequential).
#' @param pb Logical; show a progress bar for fold mapping (default \code{TRUE}).
#'
#' @details
#' For each fold:
#' \enumerate{
#' \item Subset \code{base_preds} to the analysis indices.
#' \item Fit \code{\link{fit_survmetalearner}} (time-varying NNLS) on the analysis set.
#' \item Predict stacked survival on the assessment set with
#'   \code{\link{predict_survmetalearner}}.
#' \item Compute requested metrics on the assessment set.
#' }
#'
#' Folds are \emph{stratified} by the event indicator to preserve class balance.
#' Column names in \code{base_preds} must exactly equal \code{paste0("t=", times)}.
#'
#' @section Parallelization:
#' Parallel execution is available \emph{only across folds} via \code{fmapn} when
#' \code{ncores > 1}. Base learners are assumed single-threaded; there is no inner
#' parallel layer. Choose \code{ncores} according to available CPU/RAM.
#'
#' @return An object of class \code{"cv_survmetalearner_result"} with:
#' \describe{
#'   \item{model}{Final \code{"survmetalearner"} fit on all data.}
#'   \item{cv_results}{Per-fold metric values (tibble with \code{splits}, \code{id}, \code{fold}, \code{metric}, \code{value}).}
#'   \item{summary}{Fold-aggregated mean and sd by metric.}
#'   \item{folds, metrics, ncores, pb}{Run settings for reproducibility.}
#' }
#'
#' @seealso \code{\link{fit_survmetalearner}}, \code{\link{predict_survmetalearner}},
#'   \code{\link{plot_survmetalearner_weights}}
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data(veteran, package = "survival")
#' veteran$status <- as.integer(veteran$status == 1)
#' form  <- Surv(time, status) ~ age + karno + celltype
#' times <- c(100, 200, 400)
#'
#' mod1 <- fit_coxph(form, data = veteran)
#' mod2 <- fit_rpart(form, data = veteran)
#'
#' base_models <- list(coxph = mod1, rpart = mod2)
#' base_preds  <- list(
#'   coxph = predict_coxph(mod1, veteran, times),
#'   rpart = predict_rpart(mod2, veteran, times)
#' )
#'
#' cv_res <- cv_survmetalearner(
#'   formula = form, data = veteran, times = times,
#'   base_models = base_models, base_preds = base_preds,
#'   folds = 3, metrics = c("cindex","ibs"),
#'   ncores = 2, pb = TRUE
#' )
#'
#' cv_res$summary
#' plot_survmetalearner_weights(cv_res$model)
#' }
#'
#' @export

cv_survmetalearner <- function(formula, data, times,
  base_models, base_preds,
  folds = 5,
  metrics = c("cindex", "ibs"),
  seed = 123, verbose = TRUE,
  ncores = 1, pb = TRUE) {

  stopifnot(is.list(base_models), is.list(base_preds))
  stopifnot(inherits(formula, "formula"), is.data.frame(data))
  if (!setequal(names(base_preds), names(base_models))) {
    stop("names(base_preds) must match names(base_models).")
  }

  # normalize times once
  times <- sort(unique(times))
  tnames <- paste0("t=", times)

  # attach a stable row id for indexing across rsample splits
  data <- tibble::as_tibble(data)
  data$.row_id <- seq_len(nrow(data))

  # prepare status for stratification (mirror cv_survlearner)
  tf <- terms(formula, data = data)
  outcome <- attr(tf, "variables")[[2]]
  time_col <- as.character(outcome[[2]])
  status_expr <- outcome[[3]]
  if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
    status_col <- as.character(status_expr[[2]])
    event_value <- eval(status_expr[[3]], data)
    data$.strata <- as.integer(data[[status_col]] == event_value)
  } else {
    status_col <- as.character(status_expr)
    data$.strata <- data[[status_col]]
  }

  # ensure matrices/columns align with times
  base_preds <- lapply(base_preds, function(x) {
    M <- if (!is.matrix(x)) as.matrix(x) else x
    if (!all(tnames %in% colnames(M))) {
      stop("Each base_preds matrix must have columns named exactly ", paste(tnames, collapse = ", "))
    }
    M[, tnames, drop = FALSE]
  })

  # validate row counts
  n_obs <- nrow(data)
  for (learner in names(base_preds)) {
    if (nrow(base_preds[[learner]]) != n_obs) {
      stop(sprintf("Number of rows in base_preds[['%s']] does not match data.", learner))
    }
  }

  set.seed(seed)
  
  # stratify by status like cv_survlearner
  rsplits <- rsample::vfold_cv(data, v = folds, strata = ".strata")

  # OUTER PARALLEL LAYER via fmapn
  results <- functionals::fmapn(
    list(split = rsplits$splits, id = rsplits$id, fold = seq_len(folds)),
    function(split, id, fold) {
      train <- rsample::analysis(split)
      test  <- rsample::assessment(split)

      if (verbose) cat("Processing fold", fold, "...\n")

      idx_train <- match(train$.row_id, data$.row_id)
      base_preds_train <- lapply(base_preds, function(mat) mat[idx_train, , drop = FALSE])

      # extract time/status from formula (train)
      tf <- terms(formula, data = train)
      outcome <- attr(tf, "variables")[[2]]
      time_col <- as.character(outcome[[2]])
      status_expr <- outcome[[3]]

      if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
        status_col <- as.character(status_expr[[2]])
        event_value <- eval(status_expr[[3]], train)
        status_vector <- as.integer(train[[status_col]] == event_value)
      } else {
        status_col <- as.character(status_expr)
        status_vector <- train[[status_col]]
      }

      # fit meta on training rows
      set.seed(seed + fold)
      meta_fit <- fit_survmetalearner(
        base_preds = base_preds_train,
        time = train[[time_col]],
        status = status_vector,
        times = times,
        base_models = base_models,
        formula = formula,
        data = train
      )

      # predict on test
      preds_test <- predict_survmetalearner(meta_fit, newdata = test, times = times)

      # test Surv
      time_test <- test[[time_col]]
      status_test <- if (is.call(status_expr)) {
        as.integer(test[[as.character(status_expr[[2]])]] == eval(status_expr[[3]], test))
      } else {
        test[[as.character(status_expr)]]
      }
      surv_test <- survival::Surv(time_test, status_test)

      tibble::tibble(fold = fold, metric = metrics) |>
        dplyr::mutate(value = lapply(metric, function(metric) {
          switch(metric,
            "cindex" = cindex_survmat(surv_test, predicted = preds_test, t_star = max(times)),
            "brier"  = if (length(times) != 1) stop("Brier requires single time")
                       else brier(surv_test, pre_sp = preds_test[, 1], t_star = times[[1]]),  # 🔹 scalar
            "ibs"    = ibs_survmat(surv_test, preds_test, times),
            "iae"    = iae_survmat(surv_test, preds_test, times),
            "ise"    = ise_survmat(surv_test, preds_test, times),
            stop("Unknown metric: ", metric)
          )
        })) |>
        tidyr::unnest(cols = value) |>
        dplyr::mutate(splits = list(split), id = id) |>
        dplyr::relocate(splits, id, fold)
    },
    ncores = ncores,
    pb = pb
  )

  results <- dplyr::bind_rows(results)

  # final fit on all data
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
    status_vector <- data[[status_col]]
  }

  meta_model <- fit_survmetalearner(
    base_preds = base_preds,
    time = data[[time_col]],
    status = status_vector,
    times = times,
    base_models = base_models,
    formula = formula,
    data = data
  )

  structure(list(
    model = meta_model,
    cv_results = results,
    summary = results |>
      dplyr::group_by(metric) |>
      dplyr::summarise(mean = mean(value), sd = sd(value), .groups = "drop"),
    folds = folds,
    metrics = metrics,
    ncores = ncores,
    pb = pb
  ), class = "cv_survmetalearner_result")
}
