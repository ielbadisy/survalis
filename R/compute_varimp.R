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


imp <- compute_varimp(
  model = mod_bart,
  times = c(1, 3, 5),
  metric = "ibs",
  n_repetitions = 5,
  importance_type = "delta"
)

plot_varimp(imp, use_scaled = FALSE)


plot_varimp(imp, use_scaled = TRUE)
