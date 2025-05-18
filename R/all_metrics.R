# ---- metrics-cindex.R ----
Cindex_survmat <- function(object, predicted, t_star = NULL) {
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

  C_index <- (concord + par_concord) / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}

# ---- metrics-brier.R ----
Brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  time <- object[, 1]
  status <- object[, 2]
  km_fit <- survfit(Surv(time, 1 - status) ~ 1)
  G_all <- summary(km_fit, times = sort(unique(c(t_star, time))), extend = TRUE)$surv
  names(G_all) <- as.character(sort(unique(c(t_star, time))))
  Gt <- G_all[as.character(t_star)]
  if (is.na(Gt) || Gt == 0) stop("Gt(t_star) is NA or zero")

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

  BSvalue <- mean(score_vec, na.rm = TRUE)
  names(BSvalue) <- "Brier Score"
  return(round(BSvalue, 6))
}

Brier_IBS_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")
  if (length(times) == 1) {
    return(Brier(object, pre_sp = sp_matrix[, 1], t_star = times[1]))
  }

  times <- sort(times)
  brier_vec <- vapply(seq_along(times), function(i) {
    Brier(object, pre_sp = sp_matrix[, i], t_star = times[i])
  }, numeric(1))

  dt <- diff(times)
  integral <- sum(brier_vec[-length(brier_vec)] * dt)
  IBS_value <- integral / (max(times) - min(times))
  names(IBS_value) <- "IBS"
  return(round(IBS_value, 6))
}

# ---- metrics-iaeise.R ----
IAEISE_survmat <- function(object, sp_matrix, times) {
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
  t_ISE <- (km_data$surv - pred_at_km)^2

  IAE <- sum(t_IAE[-length(t_IAE)] * dt)
  ISE <- sum(t_ISE[-length(t_ISE)] * dt)
  return(round(c(IAE = IAE, ISE = ISE), 4))
}

# ---- evaluate-survmodel.R ----
evaluate_survmodel <- function(object, sp_matrix, times,
                               metrics = c("cindex", "ibs", "iae", "ise", "brier")) {
  results <- list()
  if ("cindex" %in% metrics) {
    results$cindex <- Cindex_survmat(object, predicted = sp_matrix, t_star = max(times))
  }
  if ("brier" %in% metrics && length(times) == 1) {
    results$brier <- Brier(object, pre_sp = sp_matrix[, 1], t_star = times)
  }
  if ("ibs" %in% metrics && length(times) > 1) {
    results$ibs <- Brier_IBS_survmat(object, sp_matrix = sp_matrix, times = times)
  }
  if (any(metrics %in% c("iae", "ise"))) {
    iaeise <- IAEISE_survmat(object, sp_matrix, times)
    if ("iae" %in% metrics) results$iae <- iaeise["IAE"]
    if ("ise" %in% metrics) results$ise <- iaeise["ISE"]
  }
  return(results)
}





library(survival)
library(timereg)

# Load your functions from 'survival_metrics.R' file if needed
# source("R/survival_metrics.R")

# Example dataset
data(sTRACE)

# Step 1: Fit your functional survival model
mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)

# Step 2: Predict survival probabilities at chosen time points
time_points <- c(1, 3, 5)  # You can choose more for smoother integration
pred <- predict_aareg(mod, newdata = sTRACE, times = time_points)

# Step 3: Create Surv object for evaluation
surv_obj <- Surv(sTRACE$time, sTRACE$status == 9)

# Step 4: Evaluate all metrics
results <- evaluate_survmodel(
  object = surv_obj,
  sp_matrix = pred,
  times = time_points,
  metrics = c("cindex", "brier", "ibs", "iae", "ise")  # run all metrics
)

# Step 5: Print results
print(results)



### approach 2 

library(survival)
library(timereg)
library(purrr)
library(dplyr)
library(tibble)
library(tidyr)



# Metric names to compute
metric_list <- c("cindex", "brier", "ibs", "iae", "ise")

# Evaluate using purrr
metrics_tbl <- tibble(metric = metric_list) %>%
  mutate(value = map(metric, ~ evaluate_survmodel(
    object = surv_obj,
    sp_matrix = pred,
    times = time_points,
    metrics = .x
  )[[.x]])) %>%
  unnest(cols = value)

print(metrics_tbl)



##########



#' Evaluate survival metrics and return a tibble
#'
#' @param data A data.frame containing the survival outcome
#' @param time_col Name of column with survival time (string)
#' @param status_col Name of column with event indicator (string)
#' @param event_value Integer value indicating an event (e.g., 1 or 9)
#' @param sp_matrix Matrix/data.frame of predicted survival probabilities
#' @param times Vector of prediction time points
#' @param metrics Vector of metric names to compute (e.g. \"cindex\", \"ibs\", etc.)
#'
#' @return A tibble with columns: `metric`, `value`
#' @export
evaluate_survmodel <- function(data, time_col, status_col, event_value = 1,
  sp_matrix, times,
  metrics = c("cindex", "ibs", "iae", "ise", "brier")) {
stopifnot(is.data.frame(data))
stopifnot(time_col %in% names(data))
stopifnot(status_col %in% names(data))
stopifnot(is.matrix(sp_matrix) || is.data.frame(sp_matrix))
stopifnot(length(times) == ncol(sp_matrix))

# Convert to binary 1/0 for event/censoring
status_vector <- ifelse(data[[status_col]] == event_value, 1, 0)
surv_obj <- survival::Surv(time = data[[time_col]], event = status_vector)

tibble::tibble(metric = metrics) |>
dplyr::mutate(value = purrr::map(metric, function(metric) {
switch(metric,
"cindex" = Cindex_survmat(surv_obj, predicted = sp_matrix, t_star = max(times)),
"brier"  = {
if (length(times) != 1) stop("Brier requires a single time point")
Brier(surv_obj, pre_sp = sp_matrix[, 1], t_star = times)
},
"ibs"    = Brier_IBS_survmat(surv_obj, sp_matrix, times),
"iae"    = IAEISE_survmat(surv_obj, sp_matrix, times)["IAE"],
"ise"    = IAEISE_survmat(surv_obj, sp_matrix, times)["ISE"],
stop("Unknown metric: ", metric)
)
})) |>
tidyr::unnest(cols = value)
}


mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)
time_points <- c(1, 3, 5)
pred <- predict_aareg(mod, newdata = sTRACE, times = time_points)

evaluate_survmodel(
  data = sTRACE,
  time_col = "time",
  status_col = "status",
  event_value = 9,  # <--- HERE is the fix
  sp_matrix = pred,
  times = time_points,
  metrics = c("cindex", "ibs", "iae", "ise")
)
