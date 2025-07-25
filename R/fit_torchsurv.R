
# build DNN architecture with BatchNorm and Dropout
build_dnn <- function(input_dim, hidden, activation = "relu") {
  layers <- list()
  in_features <- input_dim
  act_fn <- switch(activation,
                   relu = nn_relu,
                   tanh = nn_tanh,
                   leaky_relu = nn_leaky_relu,
                   stop("Unsupported activation"))

  for (h in hidden) {
    layers <- append(layers, list(
      nn_linear(in_features, h),
      nn_batch_norm1d(h),
      act_fn(),
      nn_dropout(0.3)
    ))
    in_features <- h
  }
  layers <- append(layers, list(nn_linear(in_features, 1)))
  nn_sequential(!!!layers)
}

# cox partial log-likelihood loss
cox_loss <- function(pred, true) {
  time <- true[, 1]
  status <- true[, 2]
  idx <- torch_argsort(time, descending = TRUE)
  time <- time[idx]
  status <- status[idx]
  pred <- -pred[idx, 1]  # Flip sign ✅

  log_cumsum_exp <- torch_logcumsumexp(pred, dim = 1)
  event_mask <- (status == 1)
  loss <- -torch_mean(pred[event_mask] - log_cumsum_exp[event_mask])
  return(loss)
}

# internal training function
fit_coxdnn <- function(x, time, status,
                       hidden = c(32L, 16L), activation = "relu",
                       lr = 1e-4, epochs = 300L, verbose = TRUE) {
  library(torch)
  y <- cbind(time, status)
  x <- torch_tensor(as.matrix(x), dtype = torch_float())
  y <- torch_tensor(as.matrix(y), dtype = torch_float())
  net <- build_dnn(ncol(x), hidden, activation)
  optimizer <- optim_adam(net$parameters, lr = lr, weight_decay = 1e-4)

  for (epoch in 1:epochs) {
    net$train()
    optimizer$zero_grad()
    pred <- net(x)
    loss <- cox_loss(pred, y)
    loss$backward()
    optimizer$step()
    if (verbose && epoch %% 50 == 0) {
      cat("Epoch", epoch, "Loss:", loss$item(), "\n")
    }
  }

  structure(list(
    model = net, x = x, y = y,
    activation = activation, hidden = hidden,
    lr = lr, epochs = epochs
  ), class = "coxdnn")
}

# fit torchsurv model
fit_torchsurv <- function(formula, data,
                          hidden = c(32L, 32L, 32L, 16L), activation = "tanh",
                          lr = 1e-4, epochs = 1000L, verbose = TRUE) {
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(attr(mf, "terms"), data)[, -1, drop = FALSE]
  time <- y[, "time"]
  status <- y[, "status"]
  x_scaled <- scale(x)

  model <- fit_coxdnn(x_scaled, time, status, hidden, activation, lr, epochs, verbose)
  attr(x_scaled, "scaled:center") -> x_center
  attr(x_scaled, "scaled:scale") -> x_scale

  structure(list(
    model = model,
    formula = formula,
    data = data,
    xnames = colnames(x),
    x_center = x_center,
    x_scale = x_scale,
    times = sort(unique(time))
  ), class = "mlsurv_model", engine = "torchsurv")
}

# predict survival matrix
predict_torchsurv <- function(model, newdata, times) {
  library(torch)
  x <- model.matrix(delete.response(terms(model$formula)), newdata)[, model$xnames, drop = FALSE]
  x_scaled <- scale(x, center = model$x_center, scale = model$x_scale)
  x_tensor <- torch_tensor(as.matrix(x_scaled), dtype = torch_float())

  model$model$model$eval()
  with_no_grad({
    lp <- -as.numeric(model$model$model(x_tensor)[, 1])  # predicted risk
  })

  train_x <- model.matrix(delete.response(terms(model$formula)), model$data)[, model$xnames, drop = FALSE]
  train_x_scaled <- scale(train_x, center = model$x_center, scale = model$x_scale)
  train_lp <- -as.numeric(model$model$model(torch_tensor(train_x_scaled, dtype = torch_float()))[, 1])
  train_df <- data.frame(time = as.numeric(model$model$y[, 1]), status = as.numeric(model$model$y[, 2]), lp = train_lp)
  basehaz_df <- basehaz(coxph(Surv(time, status) ~ lp, data = train_df), centered = FALSE)
  H0 <- approx(basehaz_df$time, basehaz_df$hazard, xout = times, rule = 2)$y

  lp_clipped <- pmin(pmax(lp, -3), 3)
  surv_mat <- outer(lp_clipped, H0, function(lp_i, h0_j) exp(-h0_j * exp(lp_i)))
  surv_df <- as.data.frame(surv_mat)
  colnames(surv_df) <- paste0("t=", times)
  return(surv_df)
}


tune_torchsurv <- function(formula, data, times, metrics = "cindex",
                           param_grid, cv = 3, seed = 42) {
  library(purrr)
  library(dplyr)
  library(tibble)

  param_df <- tidyr::crossing(!!!param_grid)

  results <- purrr::pmap_dfr(param_df, function(hidden, lr, activation, epochs) {
    set.seed(seed)
    mod <- fit_torchsurv(formula, data, hidden = hidden,
                         lr = lr, activation = activation,
                         epochs = epochs, verbose = FALSE)

    metrics_tbl <- evaluate_survlearner(mod, metrics = metrics, times = times)
    tibble(hidden = list(hidden), lr = lr, activation = activation,
           epochs = epochs, metric = metrics_tbl$metric,
           value = metrics_tbl$value)
  })

  results %>% arrange(metric, desc(value))
}



# TEST
library(survival)
data(veteran)
mod <- fit_torchsurv(Surv(time, status) ~ age + karno + celltype, data = veteran)
pred <- predict_torchsurv(mod, newdata = veteran, times = c(30, 90, 180))
evaluate_survlearner(mod, metrics = c("cindex", "ibs", "iae", "ise"), times = c(30, 90, 180))

matplot(t(predict_torchsurv(mod, veteran[1:5, ], times = 1:365)), type = "l",
        xlab = "Time", ylab = "Survival probability", main = "Example survival curves")



### tuner
tune_torchsurv <- function(formula, data, times, metrics = "cindex",
                           param_grid, cv = 3, seed = 42) {
  library(purrr)
  library(dplyr)
  library(tibble)

  param_df <- tidyr::crossing(!!!param_grid)

  results <- purrr::pmap_dfr(param_df, function(hidden, lr, activation, epochs) {
    set.seed(seed)
    mod <- fit_torchsurv(formula, data, hidden = hidden,
                         lr = lr, activation = activation,
                         epochs = epochs, verbose = FALSE)

    metrics_tbl <- evaluate_survlearner(mod, metrics = metrics, times = times)
    tibble(hidden = list(hidden), lr = lr, activation = activation,
           epochs = epochs, metric = metrics_tbl$metric,
           value = metrics_tbl$value)
  })

  results %>% arrange(metric, desc(value))
}


grid <- list(
  hidden = list(c(16), c(32, 16)),
  lr = c(1e-4, 5e-4),
  activation = c("relu", "tanh"),
  epochs = c(300)
)

tune_torchsurv(
  Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  times = c(90),
  metrics = c("cindex", "ibs"),
  param_grid = grid
)

