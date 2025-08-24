
#' Local Surrogate Explanation for Survival Predictions (LIME-style)
#'
#' Builds a sparse, locally weighted linear surrogate model to explain a fitted
#' survival model's prediction at a specific target time for a single instance.
#' Categorical features are binarized relative to the instance of interest; local
#' weights are computed from feature-space proximity (Gower by default); and a
#' penalized (lasso) or unpenalized linear model is fit to approximate the model's
#' predicted survival probability at \code{target_time}.
#'
#' @param model An \code{mlsurv_model} fitted via \code{fit_*()}, containing a valid
#'   \code{$learner} so that \code{predict_<learner>} can be dispatched.
#' @param newdata A one-row data frame: the instance to explain (must have the same
#'   predictors as \code{baseline_data}).
#' @param baseline_data A data frame used to define the local neighborhood and to
#'   fit the surrogate (typically the model's training data).
#' @param times Numeric vector of times passed to the prediction function;
#'   must include \code{target_time} (nearest value is used).
#' @param target_time Numeric time at which to explain the prediction.
#' @param k Desired number of non-zero coefficients for the penalized surrogate
#'   (used only when \code{penalized = TRUE}); acts as a target sparsity.
#' @param dist.fun Distance function for locality weighting. Default \code{"gower"}
#'   (via \pkg{gower}); otherwise any method accepted by \code{\link[stats]{dist}}
#'   (requires \code{kernel.width}).
#' @param gower.power Power applied to \eqn{(1 - \text{gower\_dist})} when
#'   \code{dist.fun = "gower"} to control locality sharpness (default \code{5}).
#' @param kernel.width Positive numeric bandwidth used for non-Gower kernels
#'   (Gaussian weighting on pairwise distances).
#' @param penalized Logical; if \code{TRUE} (default) fits a lasso surrogate via
#'   \pkg{glmnet}; otherwise fits a weighted least squares \code{lm()}.
#' @param exclude Optional character vector of column names to exclude from the
#'   surrogate (in addition to survival outcome columns).
#'
#' @details
#' \strong{Target:} The surrogate approximates \eqn{S(t^{*} \mid x)} returned by the
#' underlying model at \code{target_time}. The nearest column in \code{times} is used.
#'
#' \strong{Weights:} Locality weights are computed from distances between
#' \code{baseline_data} rows and \code{newdata}. With \code{"gower"}, weights are
#' \eqn{(1 - \mathrm{gower\_dist})^{\mathrm{gower.power}}}; otherwise Gaussian
#' \eqn{\exp\{-d^2 / \mathrm{kernel.width}^2\}^{1/2}}.
#'
#' \strong{Sparsity:} When \code{penalized=TRUE}, a \pkg{glmnet} path is fit and the
#' solution with degrees of freedom closest to \code{k} is selected (preferring
#' exactly \code{k} if available).
#'
#' @return A data frame with one row per selected feature containing:
#' \describe{
#'   \item{feature}{Feature name (after recoding).}
#'   \item{feature_value}{Value of the feature in \code{newdata}.}
#'   \item{effect}{Local contribution \eqn{\beta_j \cdot x_j} at \code{target_time}.}
#'   \item{target_time}{The explanation time (copied for plotting).}
#' }
#' Rows are ordered by decreasing \eqn{|effect|}.
#'
#' @examples
#' \donttest{
#'  mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
#'  local_expl <- compute_surrogate(
#'  model = mod_ranger,
#'  newdata = veteran[2, ],
#'  baseline_data = veteran,
#'  times = c(100, 200, 300),
#'  target_time = 100,
#'  k = 5
#' )
#' head(local_expl)
#' }
#'
#' @export

compute_surrogate <- function(model, newdata, baseline_data,
                              times, target_time,
                              k = 5,
                              dist.fun = "gower",
                              gower.power = 5,
                              kernel.width = NULL,
                              penalized = TRUE,
                              exclude = NULL) {
  stopifnot(nrow(newdata) == 1)

  # check predict function
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model was created using fit_*() and includes a valid 'learner'.")
  }
  predict_function <- get(paste0("predict_", model$learner))

  requireNamespace("glmnet")
  requireNamespace("gower")

  # exclude time/status and user-specified vars
  response_vars <- all.vars(model$formula)[1:2]
  exclude_vars <- unique(c(response_vars, exclude))
  features <- setdiff(colnames(baseline_data), exclude_vars)

  idx_time <- which.min(abs(times - target_time))

  # predict survival probs at target time
  baseline_preds <- predict_function(model, baseline_data, times = times)
  if (is.null(dim(baseline_preds))) baseline_preds <- matrix(baseline_preds, nrow = 1)
  y <- baseline_preds[, idx_time]

  # recode X for interpretable linear modeling
  recode_data <- function(X, reference) {
    X <- as.data.frame(X)
    reference <- as.data.frame(reference)
    features <- names(reference)
    out <- as.data.frame(matrix(0, nrow = nrow(X), ncol = length(features)))
    colnames(out) <- features
    for (col in features) {
      if (is.factor(reference[[col]]) || is.character(reference[[col]])) {
        out[[col]] <- as.integer(as.character(X[[col]]) == as.character(reference[[col]]))
      } else {
        out[[col]] <- X[[col]]
      }
    }
    out
  }

  X_recode <- recode_data(baseline_data[, features, drop = FALSE], newdata[, features, drop = FALSE])
  x_interest_recode <- recode_data(newdata[, features, drop = FALSE], newdata[, features, drop = FALSE])

  # compute weights based on similarity
  # compute weights
  if (dist.fun == "gower") {
    distances <- gower_dist(baseline_data[, features, drop = FALSE],
                                   newdata[, features, drop = FALSE])
    weights <- (1 - distances)^gower.power
  } else {
    if (is.null(kernel.width)) stop("Specify kernel.width if dist.fun != 'gower'")
    dmat <- dist(rbind(newdata[, features, drop = FALSE],
                              baseline_data[, features, drop = FALSE]), method = dist.fun)
    dists <- as.matrix(dmat)[-1, 1]
    weights <- sqrt(exp(-(dists^2) / (kernel.width^2)))
  }

  # fit surrogate model (lasso or unpenalized)
  if (penalized) {
    fit <- glmnet(
      x = as.matrix(X_recode),
      y = y,
      family = "gaussian",
      weights = weights,
      standardize = TRUE,
      intercept = TRUE
    )

    df_vec <- fit$df
    if (any(df_vec == k)) {
      best_idx <- max(which(df_vec == k))
    } else {
      best_idx <- max(which(df_vec < k))
      warning("Selected fewer variables than requested k.")
    }

    coefs <- as.matrix(coef(fit, s = fit$lambda[best_idx]))
    coefs <- data.frame(feature = rownames(coefs), beta = coefs[,1])
    coefs <- coefs[coefs$feature != "(Intercept)" & coefs$beta != 0, ]

  } else {
    df <- cbind(y = y, X_recode)
    fit <- lm(y ~ ., data = df, weights = weights)
    coefs <- coef(fit)
    coefs <- data.frame(feature = names(coefs), beta = unname(coefs))
    coefs <- coefs[coefs$feature != "(Intercept)" & coefs$beta != 0, ]
  }

  # compute effects = beta * x
  effects <- sapply(coefs$feature, function(f) {
    as.numeric(x_interest_recode[[f]] * coefs$beta[coefs$feature == f])
  })

  result <- data.frame(
    feature = coefs$feature,
    feature_value = sapply(coefs$feature, function(f) as.character(newdata[[f]])),
    effect = effects,
    target_time = target_time  # added for plotting title
  )

  result <- result[order(-abs(result$effect)), ]
  rownames(result) <- NULL
  result
}


#' Plot Local Surrogate Explanation
#'
#' Visualizes the signed local effects from \code{\link{compute_surrogate}} as a
#' horizontal bar chart, optionally limiting to the top \code{n} contributors.
#'
#' @param surrogate_df A data frame returned by \code{compute_surrogate()}.
#' @param top_n Optional integer; if provided, display only the top \code{top_n}
#'   features by absolute effect.
#'
#' @return A \pkg{ggplot2} object showing feature contributions
#'   \eqn{\beta_j \cdot x_j} (positive/negative) at the target time.
#'
#' @examples
#' # plot_surrogate(local_expl, top_n = 10)
#' @export

plot_surrogate <- function(surrogate_df, top_n = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install 'ggplot2'.")
  df <- surrogate_df
  if (!is.null(top_n)) {
    df <- df[seq_len(min(top_n, nrow(df))), ]
  }

  df$feature_label <- paste0(df$feature, " = ", df$feature_value)
  df$sign <- ifelse(df$effect >= 0, "Positive", "Negative")

  target_time <- if ("target_time" %in% names(surrogate_df)) surrogate_df$target_time[1] else NA

  ggplot(df, aes(x = reorder(feature_label, abs(effect)), y = effect, fill = sign)) +
    geom_col(width = 0.6) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Local Effect (\beta * x)",
      title = paste0("Surrogate Explanation at Target Time = ", target_time)
    ) +
    scale_fill_manual(values = c("Positive" = "#4CAF50", "Negative" = "#F44336")) +
    theme_minimal() +
    theme(legend.position = "none")
}




