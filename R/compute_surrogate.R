compute_surrogate <- function(model, predict_function, newdata, baseline_data,
                              times, target_time,
                              k = 5,
                              dist.fun = "gower",
                              gower.power = 5,
                              kernel.width = NULL,
                              penalized = TRUE,
                              exclude = NULL) {
  requireNamespace("glmnet")
  requireNamespace("gower")

  stopifnot(nrow(newdata) == 1)

  # exclude time/status and user specified vars
  response_vars <- all.vars(model$formula)[1:2]
  exclude_vars <- unique(c(response_vars, exclude))
  features <- setdiff(colnames(baseline_data), exclude_vars)

  idx_time <- which.min(abs(times - target_time))

  # predict survival probs at target time
  baseline_preds <- predict_function(model, baseline_data, times = times)
  if (is.null(dim(baseline_preds))) baseline_preds <- matrix(baseline_preds, nrow = 1)
  y <- baseline_preds[, idx_time]

  # recode X: binary for factors / numeric otherwise
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

  # compute w_i based on similarity
  if (dist.fun == "gower") {
    distances <- gower::gower_dist(baseline_data[, features, drop = FALSE], newdata[, features, drop = FALSE])
    weights <- (1 - distances)^gower.power
  } else {
    if (is.null(kernel.width)) stop("Specify kernel.width if dist.fun != 'gower'")
    dists <- dist(rbind(newdata[, features, drop = FALSE], baseline_data[, features, drop = FALSE]), method = dist.fun)
    dists <- as.matrix(dists)[-1, 1]
    weights <- sqrt(exp(-(dists^2) / (kernel.width^2)))
  }

  if (penalized) {
    # lasso with penalty
    fit <- glmnet::glmnet(
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
    # unpenalized linear model
    df <- cbind(y = y, X_recode)
    fit <- lm(y ~ ., data = df, weights = weights)
    coefs <- coef(fit)
    coefs <- data.frame(feature = names(coefs), beta = unname(coefs))
    coefs <- coefs[coefs$feature != "(Intercept)" & coefs$beta != 0, ]
  }

  # compute effects = beta * X
  effects <- sapply(coefs$feature, function(f) {
    as.numeric(x_interest_recode[[f]] * coefs$beta[coefs$feature == f])
  })

  result <- data.frame(
    feature = coefs$feature,
    feature_value = sapply(coefs$feature, function(f) as.character(newdata[[f]])),
    effect = effects
  )

  result <- result[order(-abs(result$effect)), ]
  rownames(result) <- NULL
  result
}

plot_surrogate <- function(surrogate_df, top_n = NULL) {
  library(ggplot2)

  df <- surrogate_df
  if (!is.null(top_n)) {
    df <- df[seq_len(min(top_n, nrow(df))), ]
  }

  df$feature_label <- paste0(df$feature, " = ", df$feature_value)
  df$sign <- ifelse(df$effect >= 0, "Positive", "Negative")

  ggplot(df, aes(x = reorder(feature_label, abs(effect)), y = effect, fill = sign)) +
    geom_col(width = 0.6) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Local Effect (β * x)",
      title = "Surrogate Explanation at Target Time"
    ) +
    scale_fill_manual(values = c("Positive" = "#4CAF50", "Negative" = "#F44336")) +
    theme_minimal() +
    theme(legend.position = "none")
}


local_model_result <- compute_surrogate(
  model = mod_ranger,
  predict_function = predict_ranger,
  newdata = veteran[2, ],
  baseline_data = veteran,
  times = c(100, 200, 300),
  target_time = 100,
  k = 3,
  dist.fun = "gower",
  gower.power = 2,
  penalized = TRUE
)


plot_surrogate(local_model_result)
plot_surrogate(local_model_result, top_n = 3)
