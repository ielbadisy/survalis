
compute_varimp <- function(model, times,
                           metric = "ibs",
                           n_repetitions = 10,
                           seed = NULL,
                           subset = NULL) {
  if (!inherits(model, "list") || is.null(model$data)) stop("Model must include training data.")
  if (!("engine" %in% names(attributes(model)))) stop("Model must have 'engine' attribute.")

  if (!metric %in% c("cindex", "ibs", "iae", "ise", "brier")) {
    stop("Unsupported metric: ", metric)
  }

  data <- model$data
  formula <- model$formula
  pred_fun <- get(paste0("predict_", attr(model, "engine")))

  # apply subset sampling if requested
  if (!is.null(subset)) {
    if (!is.numeric(subset) || subset < 10 || subset >= nrow(data)) {
      stop("subset must be a number between 10 and nrow(data)-1")
    }
    if (!is.null(seed)) set.seed(seed)
    idx <- sample(seq_len(nrow(data)), size = subset)
    data <- data[idx, , drop = FALSE]
  }

  # extract true outcome
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

  surv_obj <- survival::Surv(time = data[[time_col]], event = status_vector)

  # base metric value
  sp_base <- pred_fun(model, newdata = data, times = times)
  base_score <- evaluate_survmodel(surv_obj, sp_base, times = times, metrics = metric)[[1]]

  # identify covariates
  exclude_vars <- c(time_col, status_col)
  covariates <- setdiff(names(data), exclude_vars)

  if (!is.null(seed)) set.seed(seed)

  results <- purrr::map(covariates, function(v) {
    scores <- replicate(n_repetitions, {
      data_perm <- data
      data_perm[[v]] <- sample(data_perm[[v]])
      sp_perm <- pred_fun(model, newdata = data_perm, times = times)
      evaluate_survmodel(surv_obj, sp_perm, times = times, metrics = metric)[[1]]
    })

    tibble::tibble(
      feature = v,
      importance = median(scores - base_score),
      importance_05 = quantile(scores - base_score, 0.05),
      importance_95 = quantile(scores - base_score, 0.95)
    )
  })

  dplyr::bind_rows(results) |>
    dplyr::arrange(desc(importance))
}


## plot varimp
plot_varimp <- function(varimp_df) {
  ggplot2::ggplot(varimp_df, ggplot2::aes(x = reorder(feature, importance), y = importance)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = importance_05, ymax = importance_95), width = 0.3) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Variable importance (diff in error)",
      title = "Permutation-based variable importance"
    ) +
    ggplot2::theme_minimal()
}

mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)

varimp_fast <- compute_varimp(
  model = mod,
  times = c(1, 3, 5),
  metric = "ibs",
  n_repetitions = 5,
  subset = 300,   # use only 300 rows to speed up | NULL for all
  seed = 123
)

plot_varimp(varimp)

