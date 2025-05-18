#' Brier Score at a single time point for survival prediction
#'
#' @param object Surv object for test set
#' @param pre_sp Vector of predicted survival probabilities at time `t_star`
#' @param t_star Time point at which to evaluate the Brier Score
#'
#' @return Brier Score (IPCW-weighted)
#' @export
Brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  time <- object[, 1]
  status <- object[, 2]

  # Estimate G(t): censoring distribution using Kaplan-Meier
  km_fit <- survfit(Surv(time, 1 - status) ~ 1)
  Gt <- summary(km_fit, times = t_star, extend = TRUE)$surv

  if (is.na(Gt)) stop("Gt(t_star) could not be estimated (too far outside range)")

  brier_sum <- 0
  n <- length(time)

  for (i in 1:n) {
    if (time[i] < t_star && status[i] == 1) {
      # Event before t_star
      Gti <- summary(km_fit, times = time[i], extend = TRUE)$surv
      if (is.na(Gti) || Gti == 0) next
      brier_sum <- brier_sum + (pre_sp[i])^2 / Gti
    } else if (time[i] >= t_star) {
      # Still at risk at t_star
      if (Gt == 0) next
      brier_sum <- brier_sum + (1 - pre_sp[i])^2 / Gt
    }
  }

  BSvalue <- brier_sum / n
  names(BSvalue) <- "Brier Score"
  return(round(BSvalue, 6))
}


#' Compute Brier Score or Integrated Brier Score (IBS) for survival prediction
#'
#' @param object A survival object from `Surv(time, status)`
#' @param sp_matrix A matrix or data.frame of predicted survival probabilities (output from `predict_aareg`)
#' @param times A single time point (for Brier Score) or vector of time points (for IBS)
#'
#' @return Brier Score or Integrated Brier Score (IBS)
#' @export
Brier_IBS_survmat <- function(object, sp_matrix, times = NULL) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (!is.matrix(sp_matrix) && !is.data.frame(sp_matrix)) stop("sp_matrix must be a matrix or data.frame")

  if (is.null(times)) {
    times <- as.numeric(gsub("t=", "", colnames(sp_matrix)))
    if (any(is.na(times))) stop("Unable to extract valid time points from column names of sp_matrix")
  }

  if (length(times) != ncol(sp_matrix)) {
    stop("Length of `times` must match number of columns in `sp_matrix`")
  }

  # Compute single Brier Score
  if (length(times) == 1) {
    pre_sp <- sp_matrix[, 1]
    t_star <- times[1]

    return(Brier(object, pre_sp = pre_sp, t_star = t_star))
  }

  # Compute IBS using trapezoidal integration
  times <- sort(times)
  brier_vec <- vapply(seq_along(times), function(i) {
    Brier(object, pre_sp = sp_matrix[, i], t_star = times[i])
  }, numeric(1))

  # Trapezoidal integration
  dt <- diff(times)
  integral <- sum(brier_vec[-length(brier_vec)] * dt)
  IBS_value <- integral / (max(times) - min(times))

  names(IBS_value) <- "IBS"
  return(round(IBS_value, 6))
}



library(survival)
data(sTRACE)

mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)

# Predict at multiple times
time_points <- c(1, 3, 5)
pred <- predict_aareg(mod, newdata = sTRACE, times = time_points)

# Brier at t = 3
Brier_IBS_survmat(Surv(sTRACE$time, sTRACE$status == 9), pred[, "t=3", drop = FALSE], times = 3)

# IBS over full range
Brier_IBS_survmat(Surv(sTRACE$time, sTRACE$status == 9), pred, times = time_points)


#**************************************

#' Vectorized Brier Score at a single time point for survival prediction
#'
#' @param object Surv object for test set
#' @param pre_sp Vector of predicted survival probabilities at time `t_star`
#' @param t_star Time point at which to evaluate the Brier Score
#'
#' @return Brier Score (IPCW-weighted)
#' @export
Brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  time <- object[, 1]
  status <- object[, 2]

  # Censoring distribution estimate
  km_fit <- survfit(Surv(time, 1 - status) ~ 1)
  G_all <- summary(km_fit, times = sort(unique(c(t_star, time))), extend = TRUE)$surv
  names(G_all) <- as.character(sort(unique(c(t_star, time))))
  Gt <- G_all[as.character(t_star)]

  if (is.na(Gt) || Gt == 0) stop("Gt(t_star) is NA or zero, can't compute Brier")

  # Vectorized computation
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



Brier_IBS_survmat(Surv(sTRACE$time, sTRACE$status == 9), pred, times = time_points)

