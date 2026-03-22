#' Permutation variable importance for survival models
#'
#' @description
#' Estimates feature importance by measuring the change in a survival metric
#' after permuting each feature.
#'
#' @param model An \code{mlsurv_model}-like object with fields \code{data}, \code{formula}, and \code{learner}.
#' @param times Numeric vector of evaluation times.
#' @param metric Character string, e.g. \code{"ibs"} or \code{"cindex"}.
#' @param n_repetitions Integer; number of permutations per feature.
#' @param seed Optional integer seed for reproducibility.
#' @param subset Optional row indices or logical vector to subset \code{model$data}.
#' @param importance_type One of \code{"delta"} (default) or \code{"mean"}; see Details.
#'
#' @return A data.frame with columns:
#' \itemize{
#'   \item \code{feature}: feature name,
#'   \item \code{value}: importance value (change in metric),
#'   \item \code{scaled_importance}: percent-scaled importance (see Details).
#' }
#'
#' @details
#' For each feature, rows are permuted \code{n_repetitions} times, predictions are recomputed,
#' and the chosen metric is compared to the baseline (unpermuted) value. The
#' \code{scaled_importance} rescales values to sum to 100%.
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#' shap_td <- compute_shap(
#'   model         = mod,
#'   newdata       = veteran[10, , drop = FALSE],
#'   baseline_data = veteran,
#'   times         = c(100),   # one time point
#'   sample.size   = 8,        # small MC
#'   aggregate     = FALSE
#' )
#' head(shap_td)
#' }
#' @export

compute_varimp <- function(model, times,
                           metric = "ibs",
                           n_repetitions = 10,
                           seed = NULL,
                           subset = NULL,
                           importance_type = c("delta", "mean")) {
  importance_type <- match.arg(importance_type)

  if (!is.list(model)) stop("Model must be a list-like object.")

  # Retrieve training data
  data <- if (!is.null(model$data)) {
    model$data
  } else if (!is.null(model$fit) && !is.null(model$fit$data)) {
    model$fit$data
  } else {
    stop("Model must include training data (e.g., model$data or model$fit$data).")
  }

  learner <- model$learner
  if (is.null(learner)) stop("Model must have a 'learner' field.")
  formula <- model$formula

  # Optional subsampling
  if (!is.null(subset)) {
    if (!is.numeric(subset) || subset < 10 || subset >= nrow(data)) {
      stop("subset must be a number between 10 and nrow(data)-1")
    }
    if (!is.null(seed)) set.seed(seed)
    data <- data[sample(seq_len(nrow(data)), subset), , drop = FALSE]
  }

  # Parse formula for time/status
  parsed_formula <- .parse_surv_formula(formula, data)
  time_col <- parsed_formula$time_col
  status_col <- parsed_formula$status_col

  if (parsed_formula$recode_status) {
    status_vector <- as.integer(data[[status_col]] == parsed_formula$event_value)
  } else {
    status_vector <- data[[status_col]]
  }

  # Compute baseline loss
  model$data <- data  # ensure correct data
  original_loss <- score_survmodel(model, times = times, metrics = metric) |>
    dplyr::filter(metric == !!metric) |>
    dplyr::pull(value)

  exclude_vars <- c(time_col, status_col)
  covariates <- setdiff(names(data), exclude_vars)
  if (!is.null(seed)) set.seed(seed)

  results <- purrr::map(covariates, function(v) {
    scores <- replicate(n_repetitions, {
      data_perm <- data
      data_perm[[v]] <- sample(data_perm[[v]])
      model$data <- data_perm  # inject permuted data
      score_survmodel(model, times = times, metrics = metric) |>
        dplyr::filter(metric == !!metric) |>
        dplyr::pull(value)
    })

    imp <- switch(importance_type,
                  "mean" = mean(scores),
                  "delta" = mean(scores - original_loss)
    )

    imp_05 <- quantile(scores, 0.05)
    imp_95 <- quantile(scores, 0.95)

    tibble::tibble(
      feature = v,
      importance = imp,
      importance_05 = imp_05,
      importance_95 = imp_95
    )
  })

  result <- dplyr::bind_rows(results)

  result$scaled_importance <- switch(importance_type,
                                     "mean"  = 100 * (result$importance - original_loss) / original_loss,
                                     "delta" = 100 * abs(result$importance) / max(abs(result$importance), na.rm = TRUE)
  )

  dplyr::arrange(result, dplyr::desc(scaled_importance))
}


#' Plot Permutation Variable Importance
#'
#' Creates a dot plot of permutation-based variable importance, using either the
#' scaled importance (default) or the raw importance column.
#'
#' @param varimp_df A data frame as returned by \code{\link{compute_varimp}}.
#' @param use_scaled Logical; if \code{TRUE} (default), plot \code{scaled_importance}
#'   (percent). If unavailable, falls back to raw \code{importance} with a warning.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' # p1 <- plot_varimp(imp, use_scaled = TRUE)
#' # p2 <- plot_varimp(imp, use_scaled = FALSE)
#' @export

plot_varimp <- function(varimp_df, use_scaled = TRUE) {
  if (use_scaled && (!"scaled_importance" %in% names(varimp_df) || all(is.na(varimp_df$scaled_importance)))) {
    warning("Scaled importance not available. Falling back to raw importance.")
    use_scaled <- FALSE
  }

  aes_x <- if (use_scaled) "scaled_importance" else "importance"

  p <- ggplot2::ggplot(varimp_df, ggplot2::aes(y = reorder(feature, !!sym(aes_x)), x = !!sym(aes_x))) +
    ggplot2::geom_point(size = 3)

  if (!use_scaled && all(c("importance_05", "importance_95") %in% names(varimp_df))) {
    # Removed as per your request: don't include error bars when using raw importance
    # p <- p + ggplot2::geom_errorbar(ggplot2::aes(xmin = importance_05, xmax = importance_95), width = 0.3)
  }

  p +
    ggplot2::labs(
      y = NULL,
      x = if (use_scaled) "Scaled importance (%)" else "Raw importance",
      title = "Permutation-based variable importance"
    ) +
    ggplot2::theme_minimal()
}
