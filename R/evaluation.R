# ---- metrics-cindex.R ----
cindex_survmat <- function(object, predicted, t_star = NULL) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")
  time <- object[, 1]
  status <- object[, 2]

  if (!is.null(t_star)) {
    t_name <- paste0("t=", t_star)
    if (!(t_name %in% colnames(predicted))) {
      stop("t_star = ", t_star, " not found in predicted survival matrix.")
    }
    surv_prob <- predicted[[t_name]]
  } else {
    surv_prob <- predicted[[ncol(predicted)]]
  }

  risk_score <- 1 - surv_prob

  permissible <- 0
  concord <- 0
  par_concord <- 0
  n <- length(time)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if ((time[i] < time[j] & status[i] == 0) | (time[j] < time[i] & status[j] == 0)) next
      if (time[i] == time[j] & status[i] == 0 & status[j] == 0) next
      permissible <- permissible + 1
      if (time[i] != time[j]) {
        if ((time[i] < time[j] & risk_score[i] > risk_score[j]) |
            (time[j] < time[i] & risk_score[j] > risk_score[i])) {
          concord <- concord + 1
        } else if (risk_score[i] == risk_score[j]) {
          par_concord <- par_concord + 0.5
        }
      } else {
        if (status[i] + status[j] > 0) {
          if (risk_score[i] == risk_score[j]) concord <- concord + 1
          else par_concord <- par_concord + 0.5
        }
      }
    }
  }

  c_index <- (concord + par_concord) / permissible
  names(c_index) <- "cindex"
  return(round(c_index, 6))
}

# ---- metrics-brier.R ----
brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  time <- object[, 1]
  status <- object[, 2]
  km_fit <- survfit(Surv(time, 1 - status) ~ 1)

  all_times <- sort(unique(c(t_star, time)))
  G_all <- summary(km_fit, times = all_times, extend = TRUE)$surv
  names(G_all) <- as.character(all_times)

  Gt <- G_all[as.character(t_star)]

  if (is.na(Gt) || Gt == 0) {
    warning("Brier: Gt(t_star) is NA or zero at t = ", t_star, "; returning NA.")
    return(NA_real_)
  }

  early_event <- which(time < t_star & status == 1)
  later_risk  <- which(time >= t_star)
  score_vec <- numeric(length(time))

  if (length(early_event) > 0) {
    Gti <- G_all[as.character(time[early_event])]
    valid <- which(!is.na(Gti) & Gti > 0)
    score_vec[early_event[valid]] <- (pre_sp[early_event[valid]]^2) / Gti[valid]
  }

  if (Gt > 0 && length(later_risk) > 0) {
    score_vec[later_risk] <- ((1 - pre_sp[later_risk])^2) / Gt
  }

  bs_value <- mean(score_vec, na.rm = TRUE)
  names(bs_value) <- "brier"
  return(round(bs_value, 6))
}

ibs_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")
  if (length(times) == 1) {
    return(brier(object, pre_sp = sp_matrix[, 1], t_star = times[1]))
  }

  times <- sort(times)
  brier_vec <- vapply(seq_along(times), function(i) {
    brier(object, pre_sp = sp_matrix[, i], t_star = times[i])
  }, numeric(1))

  dt <- diff(times)
  integral <- sum(brier_vec[-length(brier_vec)] * dt)
  ibs_value <- integral / (max(times) - min(times))
  names(ibs_value) <- "ibs"
  return(round(ibs_value, 6))
}

# ---- metrics-iaeise.R ----
iae_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  pred_at_km <- approx(x = times, y = mean_pred_surv, xout = km_data$time,
                       method = "constant", rule = 2)$y

  dt <- diff(km_data$time)
  t_IAE <- abs(km_data$surv - pred_at_km)
  iae <- sum(t_IAE[-length(t_IAE)] * dt)
  names(iae) <- "iae"
  return(round(iae, 4))
}

ise_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  pred_at_km <- approx(x = times, y = mean_pred_surv, xout = km_data$time,
                       method = "constant", rule = 2)$y

  dt <- diff(km_data$time)
  t_ISE <- (km_data$surv - pred_at_km)^2
  ise <- sum(t_ISE[-length(t_ISE)] * dt)
  names(ise) <- "ise"
  return(round(ise, 4))
}


# ---- evaluate-survmodel.R ----
evaluate_survmodel <- function(object, sp_matrix, times,
                               metrics = c("cindex", "ibs", "iae", "ise", "brier")) {
  results <- list()
  if ("cindex" %in% metrics) {
    results$cindex <- cindex_survmat(object, predicted = sp_matrix, t_star = max(times))
  }
  if ("brier" %in% metrics && length(times) == 1) {
    results$brier <- brier(object, pre_sp = sp_matrix[, 1], t_star = times)
  }
  if ("ibs" %in% metrics && length(times) > 1) {
    results$ibs <- ibs_survmat(object, sp_matrix = sp_matrix, times = times)
  }
  if ("iae" %in% metrics) {
    results$iae <- iae_survmat(object, sp_matrix, times)
  }
  if ("ise" %in% metrics) {
    results$ise <- ise_survmat(object, sp_matrix, times)
  }
  return(results)
}

evaluate_survlearner <- function(model,
                                 metrics = c("cindex", "ibs", "iae", "ise", "brier"),
                                 times) {
  if (missing(model) || !is.list(model)) stop("Must provide a fitted model object.")
  if (missing(times)) stop("Must provide time points for evaluation.")

  formula <- model$formula
  data    <- model$data
  engine  <- attr(model, "engine")

  if (is.null(engine)) stop("Model object must have an 'engine' attribute.")

  pred_fun_name <- paste0("predict_", engine)
  if (!exists(pred_fun_name, mode = "function")) {
    stop("Prediction function not found: ", pred_fun_name)
  }

  pred_fun <- get(pred_fun_name, mode = "function")
  sp_matrix <- pred_fun(model, newdata = data, times = times)

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
    event_value <- 1
    status_vector <- data[[status_col]]
  }

  surv_obj <- survival::Surv(time = data[[time_col]], event = status_vector)

  tibble::tibble(metric = metrics) |>
    dplyr::mutate(value = purrr::map(metric, function(metric) {
      switch(metric,
             "cindex" = cindex_survmat(surv_obj, predicted = sp_matrix, t_star = max(times)),
             "brier"  = {
               if (length(times) != 1) stop("Brier requires a single time point.")
               brier(surv_obj, pre_sp = sp_matrix[, 1], t_star = times)
             },
             "ibs"    = ibs_survmat(surv_obj, sp_matrix, times),
             "iae"    = iae_survmat(surv_obj, sp_matrix, times),
             "ise"    = ise_survmat(surv_obj, sp_matrix, times),
             stop("Unknown metric: ", metric)
      )
    })) |>
    tidyr::unnest(cols = value)
}

cv_survlearner <- function(formula, data,
                           fit_fun, pred_fun,
                           times,
                           metrics = c("cindex", "ibs"),
                           folds = 5,
                           seed = 123,
                           verbose = FALSE,
                           ...) {

  if ("." %in% all.vars(update(formula, . ~ 0))) {
    warning("Please avoid using '.' in the formula. Specify all predictors explicitly.")
  }

  tf <- terms(formula, data = data)
  outcome <- attr(tf, "variables")[[2]]

  if (!inherits(outcome, "call") || outcome[[1]] != as.name("Surv")) {
    stop("Formula must begin with a Surv(...) outcome.")
  }

  time_col <- as.character(outcome[[2]])
  status_expr <- outcome[[3]]

  if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
    status_col <- as.character(status_expr[[2]])
    event_value <- eval(status_expr[[3]], data)
    recode_status <- TRUE
  } else {
    status_col <- as.character(status_expr)
    event_value <- 1
    recode_status <- FALSE
  }

  all_vars <- all.vars(formula)
  n_before <- nrow(data)
  data <- tidyr::drop_na(data, dplyr::all_of(all_vars))
  n_after <- nrow(data)

  if (verbose && n_after < n_before) {
    message("Dropped ", n_before - n_after, " rows with missing values (from ", n_before, " to ", n_after, ").")
  }

  if (recode_status) {
    data[[status_col]] <- as.integer(data[[status_col]] == event_value)
  }

  set.seed(seed)
  folds_data <- rsample::vfold_cv(data, v = folds, strata = !!rlang::sym(status_col))

  purrr::pmap_dfr(
    list(split = folds_data$splits, id = folds_data$id, fold = seq_len(folds)),
    function(split, id, fold) {
      train <- rsample::analysis(split)
      test  <- rsample::assessment(split)

      model <- fit_fun(formula = formula, data = train, ...)
      pred  <- pred_fun(model, newdata = test, times = times)

      tf <- terms(formula, data = test)
      outcome <- attr(tf, "variables")[[2]]
      time_col <- as.character(outcome[[2]])
      status_expr <- outcome[[3]]

      if (is.call(status_expr) && status_expr[[1]] == as.name("==")) {
        status_col <- as.character(status_expr[[2]])
        event_value <- eval(status_expr[[3]], test)
        status_vector <- as.integer(test[[status_col]] == event_value)
      } else {
        status_col <- as.character(status_expr)
        event_value <- 1
        status_vector <- test[[status_col]]
      }

      surv_obj <- survival::Surv(time = test[[time_col]], event = status_vector)

      tibble::tibble(metric = metrics) |>
        dplyr::mutate(value = purrr::map(metric, function(metric) {
          switch(metric,
                 "cindex" = cindex_survmat(surv_obj, predicted = pred, t_star = max(times)),
                 "brier"  = {
                   if (length(times) != 1) stop("Brier requires a single time point.")
                   brier(surv_obj, pre_sp = pred[, 1], t_star = times)
                 },
                 "ibs"    = ibs_survmat(surv_obj, sp_matrix = pred, times = times),
                 "iae"    = iae_survmat(surv_obj, sp_matrix = pred, times = times),
                 "ise"    = ise_survmat(surv_obj, sp_matrix = pred, times = times),
                 stop("Unknown metric: ", metric)
          )
        })) |>
        tidyr::unnest(cols = value) |>
        dplyr::mutate(splits = list(split), id = id, fold = fold) |>
        dplyr::relocate(splits, id, fold)
    }
  )
}




cv_summary <- function(cv_results) {
  cv_results |>
    group_by(metric) |>
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd   = sd(value, na.rm = TRUE),
      n    = dplyr::n(),
      se   = sd/sqrt(n),
      lower = (mean - 1.96 * se),
      upper = (mean + 1.96 * se),
      .groups = "drop"
    )
}

cv_plot <- function(cv_results) {
  ggplot(cv_results, aes(x = metric, y = value)) +
    geom_boxplot(fill = "skyblue", outlier.shape = NA, alpha = 0.3) +
    geom_jitter(width = 0.15, alpha = 0.5) +
    labs(title = "Cross-Validation Performance", x = "Metric", y = "Value") +
    theme_minimal()
}



#----------- TEST

library(survival)
library(timereg)
library(rsample)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)


cv_results <- cv_survlearner(
  formula = Surv(time, status) ~ karno,
  data = veteran,
  fit_fun = fit_aareg,
  pred_fun = predict_aareg,
  times = c(100, 300, 500),
  metrics = c("cindex", "ibs", "iae", "ise"),
  folds = 5
)

print(cv_results)




mod <- fit_coxph(Surv(time, status) ~ age + karno, data = veteran)

evaluate_survlearner(
  model = mod,
  metrics = c("cindex", "ibs", "iae", "ise"),
  times = c(1, 5)
)

# Results
cv_summary(cv_results)
cv_plot(cv_results)







