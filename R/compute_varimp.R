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
#' @return A data.table with columns:
#' \itemize{
#'   \item \code{feature}: feature name,
#'   \item \code{importance}: importance value (change in metric),
#'   \item \code{importance_05}, \code{importance_95}: 5th/95th percentile of the
#'     per-repetition values used to compute \code{importance},
#'   \item \code{scaled_importance}: percent-scaled importance (see Details).
#' }
#' The full per-repetition values (one row per feature per repetition, columns
#' \code{feature}, \code{repetition}, \code{value}, \code{scaled_value}) are attached
#' as the \code{"raw_scores"} attribute, used by \code{\link{plot_varimp}()} to draw
#' boxplots instead of a single point per feature.
#'
#' @details
#' For each feature, rows are permuted \code{n_repetitions} times, predictions are recomputed,
#' and the chosen metric is compared to the baseline (unpermuted) value. The
#' \code{scaled_importance} rescales values to sum to 100%.
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#' imp <- compute_varimp(
#'   model = mod,
#'   times = 80,
#'   metric = "brier",
#'   n_repetitions = 3,
#'   seed = 1,
#'   subset = 40
#' )
#' head(imp)
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

  # Compute baseline loss (single metric requested, so a single row/value)
  model$data <- data  # ensure correct data
  original_loss <- score_survmodel(model, times = times, metrics = metric)$value

  exclude_vars <- c(time_col, status_col)
  covariates <- setdiff(names(data), exclude_vars)
  if (!is.null(seed)) set.seed(seed)

  per_feature <- lapply(covariates, function(v) {
    scores <- replicate(n_repetitions, {
      data_perm <- data
      data_perm[[v]] <- sample(data_perm[[v]])
      model$data <- data_perm  # inject permuted data
      score_survmodel(model, times = times, metrics = metric)$value
    })

    # values used for the point estimate are the same used for the 5/95
    # percentile band and for the raw per-repetition boxplot data, so the
    # displayed uncertainty is always about the quantity actually plotted
    values <- if (importance_type == "delta") scores - original_loss else scores

    imp <- switch(importance_type,
                  "mean" = mean(scores),
                  "delta" = mean(values)
    )

    list(
      summary = data.table::data.table(
        feature = v,
        importance = imp,
        importance_05 = stats::quantile(values, 0.05),
        importance_95 = stats::quantile(values, 0.95)
      ),
      raw = data.table::data.table(feature = v, repetition = seq_len(n_repetitions), value = values)
    )
  })

  result <- data.table::rbindlist(lapply(per_feature, `[[`, "summary"))
  raw_scores <- data.table::rbindlist(lapply(per_feature, `[[`, "raw"))

  result$scaled_importance <- switch(importance_type,
                                     "mean"  = 100 * (result$importance - original_loss) / original_loss,
                                     "delta" = 100 * abs(result$importance) / max(abs(result$importance), na.rm = TRUE)
  )
  raw_scores$scaled_value <- switch(importance_type,
                                    "mean"  = 100 * (raw_scores$value - original_loss) / original_loss,
                                    "delta" = 100 * abs(raw_scores$value) / max(abs(result$importance), na.rm = TRUE)
  )

  result <- .arrange_by_metric_dt(result, "scaled_importance", maximize = TRUE)
  attr(result, "raw_scores") <- raw_scores
  result
}


#' Plot Permutation Variable Importance
#'
#' Plots the distribution of per-repetition permutation importance values as a
#' boxplot per feature, using either the scaled importance (default) or the raw
#' importance column. Falls back to a single point per feature (with no
#' distribution shown) for \code{varimp_df} objects that lack the
#' \code{"raw_scores"} attribute (e.g., hand-built summary tables).
#'
#' @param varimp_df A data frame as returned by \code{\link{compute_varimp}}.
#' @param use_scaled Logical; if \code{TRUE} (default), plot \code{scaled_importance}
#'   (percent). If unavailable, falls back to raw \code{importance} with a warning.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely (e.g., for journals
#'   requiring caption-only figures).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#' imp <- compute_varimp(
#'   model = mod,
#'   times = 80,
#'   metric = "brier",
#'   n_repetitions = 3,
#'   seed = 1,
#'   subset = 40
#' )
#' plot_varimp(imp, use_scaled = TRUE)
#' plot_varimp(imp, use_scaled = FALSE)
#' @export

plot_varimp <- function(varimp_df, use_scaled = TRUE, title) {
  if (use_scaled && (!"scaled_importance" %in% names(varimp_df) || all(is.na(varimp_df$scaled_importance)))) {
    warning("Scaled importance not available. Falling back to raw importance.")
    use_scaled <- FALSE
  }
  if (missing(title)) title <- "Permutation-based variable importance"

  aes_x <- if (use_scaled) "scaled_importance" else "importance"
  feature_levels <- varimp_df$feature[order(varimp_df[[aes_x]])]
  xlab <- if (use_scaled) "Scaled importance (%)" else "Raw importance"

  raw <- attr(varimp_df, "raw_scores")

  if (!is.null(raw)) {
    plot_df <- data.table::copy(raw)
    plot_df$feature <- factor(plot_df$feature, levels = feature_levels)
    value_col <- if (use_scaled) "scaled_value" else "value"

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(y = feature, x = .data[[value_col]])) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
        ggplot2::geom_boxplot(fill = .survalis_palette[1], alpha = 0.6, outlier.alpha = 0.4) +
        ggplot2::labs(y = NULL, x = xlab) +
        theme_survalis()
    if (!is.null(title)) p <- p + ggplot2::labs(title = title)
    return(p)
  }

  p <- ggplot2::ggplot(varimp_df, ggplot2::aes(y = reorder(feature, !!sym(aes_x)), x = !!sym(aes_x))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_point(size = 3, color = .survalis_palette[1]) +
    ggplot2::labs(y = NULL, x = xlab) +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}
