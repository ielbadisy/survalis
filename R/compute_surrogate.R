
## NEED TO MAKE THESE ADJUSTMENTS

### penalization should be an option (defaut?)
### make sure to exclude time and status
compute_surrogate <- function(model, predict_function, newdata, baseline_data,
                                 times, target_time,
                                 k = 5,
                                 dist.fun = "gower",
                                 gower.power = 5,
                                 kernel.width = NULL) {

  ## fix this: avoid heavy loading, need to do it through the DECRIPTION
  requireNamespace("glmnet")
  requireNamespace("gower")
  requireNamespace("data.table")

  stopifnot(nrow(newdata) == 1)

  #we exclude time and status from features
  response_vars <- all.vars(model$formula)[1:2]  # time and status
  features <- setdiff(colnames(baseline_data), response_vars)

  idx_time <- which.min(abs(times - target_time))

  # make predictions
  baseline_preds <- predict_function(model, baseline_data, times = times)
  if (is.null(dim(baseline_preds))) baseline_preds <- matrix(baseline_preds, nrow = 1)
  y <- baseline_preds[, idx_time]

  # recoding
  recode_data <- function(X, reference) {
    X <- as.data.frame(X)
    reference <- as.data.frame(reference)
    features <- names(reference)
    out <- matrix(0, nrow = nrow(X), ncol = length(features))
    colnames(out) <- features
    out <- as.data.frame(out)

    for (col in features) {
      if (is.factor(reference[[col]]) || is.character(reference[[col]])) {
        out[[col]] <- as.integer(as.character(X[[col]]) == as.character(reference[[col]]))
      } else {
        out[[col]] <- X[[col]]
      }
    }

    out
  }

  X_recode <- recode_data(baseline_data, newdata)
  x_interest_recode <- recode_data(newdata, newdata)

  # compute proximity weights (default using gower dist)
  if (dist.fun == "gower") {
    distances <- gower::gower_dist(baseline_data[, features, drop = FALSE], newdata[, features, drop = FALSE])
    weights <- (1 - distances)^gower.power
  } else {
    if (is.null(kernel.width)) {
      stop("You must specify kernel.width if dist.fun is not 'gower'")
    }
    dists <- dist(rbind(newdata[, features, drop = FALSE], baseline_data[, features, drop = FALSE]), method = dist.fun)
    dists <- as.matrix(dists)[-1, 1]
    weights <- sqrt(exp(-(dists^2) / (kernel.width^2)))
  }

  # fit weighted Lasso
  fit <- glmnet::glmnet(
    x = as.matrix(X_recode),
    y = y,
    family = "gaussian",
    weights = weights,
    standardize = TRUE,
    intercept = TRUE
  )

  # select best lambda -> penalize for more clair interpretations (need to justfy this choice)
  df_vec <- fit$df
  if (any(df_vec == k)) {
    best_idx <- max(which(df_vec == k))
  } else {
    best_idx <- max(which(df_vec < k))
    warning("Selected fewer variables than requested k")
  }

  coefs <- as.matrix(coef(fit, s = fit$lambda[best_idx]))
  coefs <- data.frame(feature = rownames(coefs), beta = coefs[,1])
  coefs <- coefs[coefs$feature != "(Intercept)", ]
  coefs <- coefs[coefs$beta != 0, ]

  # calculate effect (beta * feature value)
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



local_model_result <- local_surrogate_surv(
  model = mod_ranger,
  predict_function = predict_ranger,
  newdata = veteran[1, ],
  baseline_data = veteran,
  times = times,
  target_time = 100,
  k = 3,
  dist.fun = "gower",
  gower.power = 2
)

print(local_model_result)


#----- plot_surrogate()
plot_surrogate <- function(local_model_result, top_n = NULL) {
  stopifnot(is.data.frame(local_model_result))

  plot_data <- local_model_result
  if (!is.null(top_n)) {
    top_features <- head(plot_data[order(-abs(plot_data$effect)), ], n = top_n)
    plot_data <- top_features
  }

  plot_data$feature_label <- paste0(plot_data$feature, " = ", plot_data$feature_value)

  p <- ggplot(plot_data, aes(x = reorder(feature_label, effect), y = effect)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
      title = "Local Surrogate Model: Feature Effects",
      x = "Feature (value)",
      y = "Effect on Survival Probability"
    ) +
    theme_minimal(base_size = 14)

  return(p)
}

library(ggplot2)
plot_surrogate(local_model_result)


plot_surrogate(local_model_result, top_n = 2)
