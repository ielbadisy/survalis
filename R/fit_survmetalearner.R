
#' Fit a Stacked Survival Meta‑Learner (Time‑Varying NNLS)
#'
#' Learns nonnegative, time‑specific stacking weights over a set of base survival
#' learners by solving a nonnegative least squares (NNLS) problem at each
#' requested time \eqn{t^{*}}. For each \eqn{t^{*}}, the target is the indicator
#' \eqn{Y = I(T > t^{*})}; the features are the base learners' predicted survival
#' probabilities \eqn{S_\ell(t^{*} \mid x)}. Weights are constrained to be
#' nonnegative and are normalized to sum to 1 per time point.
#'
#' @param base_preds A named list of matrices/data frames, one per base learner,
#'   each of dimension \code{n x length(times)}, with columns named \code{"t=<time>"}.
#' @param time Numeric vector of observed event/censoring times (length \code{n}).
#' @param status Numeric/binary vector of event indicators (1=event, 0=censor) (length \code{n}).
#' @param times Numeric vector of evaluation times at which to learn weights.
#' @param base_models A named list of fitted base learner objects; names must
#'   match \code{names(base_preds)} and the learner names used by \code{predict_*()}.
#' @param formula A \code{Surv()} formula that was used for the base models
#'   (stored for metadata and downstream scoring).
#' @param data The training data frame used for the base models (stored for metadata).
#'
#' @details
#' For each \code{t} in \code{times}, this function fits
#' \deqn{\min_{w \ge 0} \\mid Y - S w \|_2^2 \quad \text{s.t.} \ \sum_j w_j = 1}
#' where \eqn{Y = I(T > t)} and \eqn{S} is the matrix of base survival
#' probabilities at \eqn{t}. The NNLS solution from \pkg{nnls} is renormalized
#' to sum to 1 (if the solution is all zeros, weights remain \code{NA} for that time).
#'
#' @return An object of class \code{c("mlsurv_model","survmetalearner")} with elements:
#' \describe{
#'   \item{weights}{Matrix \code{L x T} of nonnegative stacking weights, rows=learners, cols=\code{"t=<time>"}.}
#'   \item{base_models}{Named list of base learner fits.}
#'   \item{base_preds}{Named list of base prediction matrices on training data.}
#'   \item{learners}{Character vector of learner names (from \code{names(base_preds)}).}
#'   \item{formula, data, time, status}{Training metadata for scoring/reporting.}
#'   \item{learner}{The string \code{"survmetalearner"} (for \code{predict_*} dispatch).}
#' }
#'
#' @seealso \code{\link{predict_survmetalearner}}, \code{\link{plot_survmetalearner_weights}},
#'   \code{\link{cv_survmetalearner}}
#'
#' @examples
#' # See cv_survmetalearner() and predict_survmetalearner() examples.
#' @export

fit_survmetalearner <- function(base_preds, time, status, times,
                            base_models, formula, data) {
  stopifnot(is.list(base_preds), is.list(base_models), is.data.frame(data))

  base_preds <- lapply(base_preds, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })

  learners <- names(base_preds)
  T <- length(times)
  weight_matrix <- matrix(NA_real_, nrow = length(learners), ncol = T,
                          dimnames = list(learners, paste0("t=", times)))

  for (j in seq_along(times)) {
    t_star <- times[j]
    Y <- as.integer(time > t_star)
    X <- sapply(base_preds, function(S) S[, j])
    fit <- nnls::nnls(as.matrix(X), Y)
    w <- coef(fit); w <- w / sum(w)
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
    learner = "survmetalearner",         # for predict_* lookup
    learners = learners
  ), class = c("mlsurv_model", "survmetalearner"), engine = "survalis")  # dual class

}

#' Predict with a Stacked Survival Meta‑Learner
#'
#' Produces stacked survival probabilities by combining base learner predictions
#' via the time‑specific weights learned by \code{\link{fit_survmetalearner}}.
#'
#' @param model A \code{"survmetalearner"} object returned by \code{fit_survmetalearner()}.
#' @param newdata A data frame of new observations for prediction.
#' @param times Numeric vector of evaluation times (must be a subset of the times
#'   used to train the meta‑learner).
#'
#' @details
#' For each base learner listed in \code{model$learners}, the corresponding
#' \code{predict_<learner>} function is called to obtain \eqn{S_\ell(t \mid x)}.
#' Stacked predictions are computed as \eqn{\hat{S}(t \mid x) = \sum_\ell w_\ell(t) S_\ell(t \mid x)},
#' where \eqn{w_\ell(t)} are the learned nonnegative weights for time \eqn{t}.
#'
#' @return A data frame with one row per observation and one column per requested
#'   time (columns named \code{"t=<time>"}), containing stacked survival probabilities.
#'
#' @seealso \code{\link{fit_survmetalearner}}, \code{\link{plot_survmetalearner_weights}}
#'
#' @examples
#' # pred <- predict_survmetalearner(meta_model, newdata = df_test, times = c(200, 500, 800))
#' @export

predict_survmetalearner <- function(model, newdata, times) {
  stopifnot(inherits(model, "survmetalearner"))

  learners <- model$learners
  base_models <- model$base_models
  requested_tnames <- paste0("t=", times)
  W <- model$weights[, requested_tnames, drop = FALSE]

  base_preds_test <- lapply(learners, function(learner) {
    pred_fun <- get(paste0("predict_", learner), mode = "function")
    pred <- pred_fun(base_models[[learner]], newdata = newdata, times = times)
    pred <- as.matrix(pred)
    colnames(pred) <- requested_tnames
    pred
  })
  names(base_preds_test) <- learners

  n <- nrow(base_preds_test[[1]])
  T <- length(times)
  pred_array <- array(NA_real_, dim = c(n, length(learners), T))
  for (k in seq_along(learners)) {
    pred_array[, k, ] <- base_preds_test[[k]]
  }

  survmat <- matrix(NA_real_, nrow = n, ncol = T)
  for (j in seq_len(T)) {
    survmat[, j] <- pred_array[, , j] %*% W[, j]
  }

  .finalize_survmat(survmat, times = times)
}

#' Plot Time‑Varying Stacking Weights
#'
#' Visualizes the learned nonnegative NNLS stacking weights \eqn{w_\ell(t)} over time
#' for each base learner.
#'
#' @param model A \code{"survmetalearner"} object from \code{\link{fit_survmetalearner}}.
#'
#' @return A \pkg{ggplot2} object showing weight trajectories (one line per learner).
#'
#' @examples
#' # plot_survmetalearner_weights(meta_model)
#' @export

plot_survmetalearner_weights <- function(model) {
  stopifnot(inherits(model, "survmetalearner"))

  W <- as.data.frame(model$weights)
  W$learner <- rownames(model$weights)

  W_long <- W |>
    pivot_longer(-learner, names_to = "time", values_to = "weight") |>
    mutate(time = as.numeric(sub("t=", "", time)))

  ggplot(W_long, aes(x = time, y = weight, color = learner)) +
    geom_line(linewidth = 1.2) +
    labs(x = "Time", y = "Weight", title = "NNLS stacking weights over time") +
    theme_minimal(base_size = 14)
}

#' Cross‑Validate a Stacked Survival Meta‑Learner
#'
#' Performs v‑fold cross‑validation for the NNLS stacking meta‑learner over
#' a fixed set of base learners and their predictions. On each fold, stacking
#' weights are learned on the analysis set and evaluated on the assessment set.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing all variables referenced in \code{formula}.
#' @param times Numeric vector of evaluation times for stacking and scoring.
#' @param base_models A named list of fitted base learner models (used to predict
#'   on assessment folds and for the final refit on full data).
#' @param base_preds A named list of \emph{training‑set} prediction matrices
#'   (rows align with \code{data}, columns align with \code{times}); names must
#'   match \code{base_models}.
#' @param folds Integer; number of CV folds (default \code{5}).
#' @param metrics Character vector of metrics to compute (default \code{c("cindex","ibs")}).
#'   Supported: \code{"cindex"}, \code{"brier"} (single time), \code{"ibs"},
#'   \code{"iae"}, \code{"ise"}, \code{"ece"} (single time).
#' @param seed Integer random seed (default \code{123}).
#' @param verbose Logical; if \code{TRUE}, prints fold progress (default \code{TRUE}).
#'
#' @details
#' For each fold: (1) subset \code{base_preds} to the analysis indices;
#' (2) learn time‑specific NNLS weights with \code{\link{fit_survmetalearner}};
#'
#' (3) predict stacked survival on the assessment set with \code{\link{predict_survmetalearner}};
#' (4) compute requested metrics. After CV, a final meta‑learner is fit on the full data.
#'
#' @return An object of class \code{"cv_survmetalearner_result"} with components:
#' \describe{
#'   \item{model}{Final \code{"survmetalearner"} fit on all data.}
#'   \item{cv_results}{Per‑fold metric values (tibble).}
#'   \item{summary}{Fold‑aggregated mean and sd by metric (tibble).}
#'   \item{folds, metrics}{CV settings.}
#' }
#'
#' @seealso \code{\link{fit_survmetalearner}}, \code{\link{predict_survmetalearner}},
#'   \code{\link{plot_survmetalearner_weights}}
#'
#' @examples
#' \dontrun{
#' form  <- Surv(time, status) ~ age + karno + trt + diagtime + celltype
#' times <- c(200, 500, 800)
#'
#' # Fit a few base learners (examples)
#' mod_rpart   <- fit_rpart(form, data = veteran)
#' mod_ranger  <- fit_ranger(form, data = veteran)
#' mod_cforest <- fit_cforest(form, data = veteran)
#'
#' base_models <- list(rpart = mod_rpart, ranger = mod_ranger, cforest = mod_cforest)
#' base_preds  <- list(
#'   rpart   = predict_rpart(mod_rpart,   veteran, times),
#'   ranger  = predict_ranger(mod_ranger, veteran, times),
#'   cforest = predict_cforest(mod_cforest, veteran, times)
#' )
#'
#' cv_res <- cv_survmetalearner(
#'   formula = form, data = veteran, times = times,
#'   base_models = base_models, base_preds = base_preds,
#'   folds = 3, metrics = c("cindex","ibs")
#' )
#'
#' cv_res$summary
#' plot_survmetalearner_weights(cv_res$model)
#' }
#' @export

cv_survmetalearner <- function(formula, data, times,
                           base_models, base_preds,
                           folds = 5,
                           metrics = c("cindex", "ibs"),
                           seed = 123, verbose = TRUE) {
  stopifnot(is.list(base_models), is.list(base_preds))
  stopifnot(inherits(formula, "formula"), is.data.frame(data))

  # Convert all base_preds to matrices (if needed)
  base_preds <- lapply(base_preds, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })

  # Validate shape
  n_obs <- nrow(data)
  for (learner in names(base_preds)) {
    if (nrow(base_preds[[learner]]) != n_obs) {
      stop(sprintf("Number of rows in base_preds[['%s']] does not match data.", learner))
    }
  }

  set.seed(seed)
  rsplits <- rsample::vfold_cv(data, v = folds)

  results <- purrr::map_dfr(seq_along(rsplits$splits), function(i) {
    split <- rsplits$splits[[i]]
    train <- rsample::analysis(split)
    test  <- rsample::assessment(split)

    if (verbose) cat("Processing fold", i, "...\n")

    # subset base_preds for training rows
    idx_train <- as.integer(rownames(data)) %in% as.integer(rownames(train))
    base_preds_train <- lapply(base_preds, function(mat) mat[idx_train, , drop = FALSE])

    # extract time/status from formula
    parsed_formula <- .parse_surv_formula(formula, train)
    time_col <- parsed_formula$time_col
    status_col <- parsed_formula$status_col

    if (parsed_formula$recode_status) {
      status_vector <- as.integer(train[[status_col]] == parsed_formula$event_value)
    } else {
      status_vector <- train[[status_col]]
    }

    # fit survmetalearner on training fold
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

    # extract test survival object
    time_test <- test[[time_col]]
    status_test <- if (parsed_formula$recode_status) {
      as.integer(test[[status_col]] == parsed_formula$event_value)
    } else {
      test[[status_col]]
    }
    surv_test <- survival::Surv(time_test, status_test)

    # score metrics
    tibble::tibble(
      fold = i,
      metric = metrics
    ) |>
      dplyr::mutate(value = purrr::map(metric, function(metric) {
        switch(metric,
               "cindex" = cindex_survmat(surv_test, predicted = preds_test, t_star = max(times)),
               "auc"    = auc_survmat(surv_test, predicted = preds_test, t_star = max(times)),
               "brier"  = if (length(times) != 1) stop("Brier requires single time") else
                 brier(surv_test, pre_sp = preds_test[, 1], t_star = times),
               "ibs"    = ibs_survmat(surv_test, preds_test, times),
               "iae"    = iae_survmat(surv_test, preds_test, times),
               "ise"    = ise_survmat(surv_test, preds_test, times),
               "ece"    = if (length(times) != 1) stop("ECE requires single time") else
                 ece_survmat(surv_test, preds_test, t_star = times),
               stop("Unknown metric: ", metric)
        )
      })) |>
      tidyr::unnest(cols = value)
  })

  # Fit final survmetalearner on full data
  parsed_formula <- .parse_surv_formula(formula, data)
  time_col <- parsed_formula$time_col
  status_col <- parsed_formula$status_col

  if (parsed_formula$recode_status) {
    status_vector <- as.integer(data[[status_col]] == parsed_formula$event_value)
  } else {
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
    summary = results |> dplyr::group_by(metric) |> dplyr::summarise(mean = mean(value), sd = sd(value)),
    folds = folds,
    metrics = metrics
  ), class = "cv_survmetalearner_result")
}
