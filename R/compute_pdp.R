#' Partial Dependence and ICE for Survival Predictions
#'
#' Computes partial dependence (PDP) and/or individual conditional expectation (ICE)
#' curves of predicted survival probabilities for a single feature at one or more
#' evaluation times. Works with any learner fitted via `fit_*()` that exposes a
#' matching `predict_*()` method returning survival probabilities.
#'
#' @param model An `mlsurv_model` created by a `fit_*()` function; must contain a
#'   valid `learner` so the appropriate `predict_<learner>()` can be dispatched.
#' @param data A data frame used to construct PDP/ICE profiles (typically the
#'   training data).
#' @param feature Character scalar; the feature name to analyze (numeric or categorical).
#' @param times Numeric vector of evaluation times at which survival probabilities
#'   are computed.
#' @param method One of `"pdp"`, `"ice"`, or `"pdp+ice"` (default). Controls which
#'   profiles are produced.
#' @param grid.size Integer number of grid points for numeric features (default 20).
#'   Ignored for categorical features (levels are used).
#'
#' @details
#' For numeric features, a regular grid over the observed range is used; for
#' categorical features, all observed levels are used. For each grid value,
#' predictions are made for every row of `data` with the feature forced to the
#' grid value, yielding ICE curves per row and PDP as the average across rows.
#'
#' If multiple `times` are supplied, outputs are stacked in long format with a
#' `time` column; an additional integrated PDP is computed via trapezoidal rule
#' when more than one time is provided.
#'
#' @return A list with elements:
#' \describe{
#'   \item{results}{Long data frame with columns:
#'     \code{surv_prob}, \code{time}, \code{type} (pdp/ice),
#'     \code{.id} (row id for ICE), and the analyzed \code{feature}.}
#'   \item{pdp_integrated}{(Optional) Data frame with \code{feature} and
#'     \code{integrated_surv} (time-integrated PDP), present only if
#'     \code{length(times) > 1} and PDP was requested.}
#' }
#'
#' @examples
#' # mod_ranger <- fit_ranger(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' # pdp_age <- compute_pdp(mod_ranger, data = veteran, feature = "age",
#' #                        times = c(100, 200, 300), method = "pdp+ice")
#' # str(pdp_age)
#'
#' @export

compute_pdp <- function(model, data, feature, times,
                        method = "pdp+ice", grid.size = 20) {
  # Automatically determine prediction function
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model has a valid `learner` and was created using `fit_*()`.")
  }
  predict_function <- get(paste0("predict_", model$learner))

  if (!feature %in% names(data)) stop("Feature not found in data.")

  is_categorical <- is.factor(data[[feature]]) || is.character(data[[feature]])

  if (is_categorical) {
    grid_values <- sort(unique(data[[feature]]))
  } else {
    feature_range <- range(data[[feature]], na.rm = TRUE)
    grid_values <- seq(feature_range[1], feature_range[2], length.out = grid.size)
  }

  n_obs <- nrow(data)
  n_grid <- length(grid_values)
  expanded_data <- data[rep(seq_len(n_obs), each = n_grid), , drop = FALSE]
  expanded_data[[feature]] <- rep(grid_values, times = n_obs)

  surv_preds <- predict_function(model, newdata = expanded_data, times = times)
  if (is.null(dim(surv_preds))) surv_preds <- matrix(surv_preds, ncol = 1)
  if (is.null(colnames(surv_preds))) colnames(surv_preds) <- paste0("t=", times)

  surv_df <- as.data.frame(surv_preds)
  colnames(surv_df) <- paste0("t=", times)
  surv_df[[feature]] <- rep(grid_values, times = n_obs)
  surv_df$.id <- rep(seq_len(n_obs), each = n_grid)

  results_list <- list()

  if (method %in% c("ice", "pdp+ice")) {
    ice_long <- reshape(surv_df,
                        varying = paste0("t=", times),
                        v.names = "surv_prob",
                        timevar = "time",
                        times = times,
                        idvar = c(".id", feature),
                        direction = "long")
    ice_long$type <- "ice"
    rownames(ice_long) <- NULL
    results_list$ice <- ice_long
  }

  if (method %in% c("pdp", "pdp+ice")) {
    means <- aggregate(surv_df[, paste0("t=", times), drop = FALSE],
                       by = list(surv_df[[feature]]),
                       FUN = mean)
    colnames(means)[1] <- feature
    pdp_long <- reshape(means,
                        varying = paste0("t=", times),
                        v.names = "surv_prob",
                        timevar = "time",
                        times = times,
                        idvar = feature,
                        direction = "long")
    pdp_long$type <- "pdp"
    pdp_long$.id <- NA
    rownames(pdp_long) <- NULL
    results_list$pdp <- pdp_long
  }

  final_result <- do.call(rbind, results_list)
  final_result$time <- as.numeric(as.character(final_result$time))

  pdp_integrated <- NULL
  if ("pdp" %in% final_result$type && length(unique(final_result$time)) > 1) {
    integrate_fun <- function(x, y) {
      x <- as.numeric(x)
      y <- as.numeric(y)
      sum((head(y, -1) + tail(y, -1)) / 2 * diff(x))
    }

    split_pdp <- split(final_result[final_result$type == "pdp", ], final_result[final_result$type == "pdp", feature])
    pdp_integrated <- do.call(rbind, lapply(split_pdp, function(df) {
      int_surv <- integrate_fun(df$time, df$surv_prob)
      f_val <- unique(df[[feature]])
      out <- data.frame(integrated_surv = int_surv, stringsAsFactors = FALSE)
      out[[feature]] <- f_val
      out
    }))
    rownames(pdp_integrated) <- NULL
    pdp_integrated <- pdp_integrated[, c(feature, "integrated_surv")]
  }

  list(
    results = final_result,
    pdp_integrated = pdp_integrated
  )
}



#' Plot PDP/ICE Curves for Survival Models
#'
#' Plots partial dependence (PDP) and/or individual conditional expectation (ICE)
#' results returned by \code{\link{compute_pdp}} either per evaluation time or as
#' an integrated PDP over time.
#'
#' @param pdp_ice_output The list returned by \code{compute_pdp()}.
#' @param feature Character scalar; the same feature analyzed in \code{compute_pdp()}.
#' @param method One of `"pdp"`, `"ice"`, or `"pdp+ice"` indicating which curves
#'   to draw in per-time plots (ignored for integrated plots).
#' @param ids Optional vector of row ids (`.id`) to subset ICE curves for clarity
#'   in per-time plots.
#' @param which One of `"per_time"` (facet by time) or `"integrated"` (plot the
#'   time-integrated PDP if available).
#' @param alpha_ice Alpha transparency for ICE lines/boxes in per-time plots (default 0.2).
#' @param smooth Logical; if `TRUE` and the feature is numeric, apply a smooth
#'   curve (`loess`) for integrated PDP; otherwise draw a line.
#'
#' @details
#' \strong{Per-time:} For numeric features, draws ICE lines and PDP overlays per time.
#' For categorical features, shows ICE as boxplots per level and PDP as point summaries.
#' \strong{Integrated:} Plots the PDP integrated across time (if provided by
#' \code{compute_pdp()}); numeric features can be smoothed with \code{smooth=TRUE}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' # plot_pdp(pdp_age, feature = "age", which = "per_time")
#' # plot_pdp(pdp_age, feature = "age", which = "integrated", smooth = TRUE)
#'
#' @export

plot_pdp <- function(pdp_ice_output, feature,
                     method = "pdp+ice", ids = NULL,
                     which = c("per_time", "integrated"),
                     alpha_ice = 0.2, smooth = FALSE) {
  requireNamespace("ggplot2")
  requireNamespace("data.table")

  which <- match.arg(which)
  results <- as.data.table(pdp_ice_output$results)
  feature_data <- results[[feature]]
  is_categorical <- is.factor(feature_data) || is.character(feature_data)

  if (which == "integrated") {

    if (is.null(pdp_ice_output$pdp_integrated)) stop("No integrated PDP found.")

    if (!"integrated_surv" %in% names(pdp_ice_output$pdp_integrated)) stop("Invalid integrated PDP.")

    ice_data <- results[type == "ice"]
    ice_integrated <- ice_data[, .(
      integrated_surv = sum((head(surv_prob, -1) + tail(surv_prob, -1)) / 2 * diff(head(time, -1))),
      time_n = .N
    ), by = c(".id", feature)]

    pdp_integrated <- pdp_ice_output$pdp_integrated

    if (is_categorical) {
      p <- ggplot() +
        geom_boxplot(data = ice_integrated, aes(x = .data[[feature]], y = integrated_surv),
                              alpha = alpha_ice, fill = "pink") +
        geom_point(data = pdp_integrated, aes(x = .data[[feature]], y = integrated_surv),
                            shape = 21, size = 3, fill = "black") +
        theme_minimal(base_size = 13) +
        labs(
          title = paste("Integrated PDP with ICE Boxplot for", feature),
          x = feature,
          y = "Integrated Survival"
        )
      return(p)
    } else {
      p <- ggplot(pdp_integrated, aes(x = .data[[feature]], y = integrated_surv)) +
        theme_minimal(base_size = 13) +
        labs(
          title = paste("Integrated PDP over Survival Time:", feature),
          x = feature,
          y = "Integrated Survival"
        )
      if (smooth) {
        p <- p + geom_smooth(method = "loess", se = FALSE, color = "steelblue", linewidth = 1.2)
      } else {
        p <- p + geom_line(color = "steelblue", linewidth = 1.2)
      }
      return(p)
    }
  }

  # per-time PDP/ICE plot
  plot_data <- results[type %in% switch(method,
                                        "pdp" = "pdp",
                                        "ice" = "ice",
                                        "pdp+ice" = c("pdp", "ice"),
                                        stop("Invalid method"))]

  if (!is.null(ids) && "ice" %in% plot_data$type) {
    plot_data <- plot_data[type != "ice" | .id %in% ids]
  }

  if (is_categorical) {
    p <- ggplot(plot_data, aes(x = .data[[feature]], y = surv_prob, fill = type)) +
      theme_minimal(base_size = 13) +
      facet_wrap(~ time, scales = "free_y") +
      labs(
        title = paste("Survival", toupper(method), "for", feature),
        x = feature,
        y = "Survival Probability",
        fill = "Type"
      )

    if ("ice" %in% plot_data$type) {
      p <- p + geom_boxplot(data = plot_data[type == "ice"], alpha = alpha_ice, position = "dodge")
    }

    if ("pdp" %in% plot_data$type) {
      p <- p + stat_summary(data = plot_data[type == "pdp"],
                                     fun = mean, geom = "point",
                                     shape = 21, fill = "black", size = 3,
                                     position = position_dodge(width = 0.75))
    }

  } else {
    p <- ggplot(
      plot_data,
      aes(
        x = .data[[feature]],
        y = surv_prob,
        group = interaction(.id, time),
        color = type
      )
    ) +
      theme_minimal(base_size = 13) +
      facet_wrap(~ time, scales = "free_y") +
      labs(
        title = paste("Survival", toupper(method), "for", feature),
        x = feature,
        y = "Survival Probability",
        color = NULL
      ) +
      guides(color = guide_legend(override.aes = list(alpha = 1)))

    if ("ice" %in% plot_data$type) {
      p <- p + geom_line(data = plot_data[type == "ice"], alpha = alpha_ice)
    }

    if ("pdp" %in% plot_data$type) {
      p <- p + geom_line(data = plot_data[type == "pdp"], linewidth = 1.2)
    }
  }

  return(p)
}



