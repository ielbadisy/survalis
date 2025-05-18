#' Integrated Absolute Error (IAE) and Integrated Squared Error (ISE)
#' for survival prediction based on predicted survival probabilities
#'
#' @param object A survival object from `Surv(time, status)`
#' @param sp_matrix A matrix or data.frame of predicted survival probabilities (rows = obs, cols = times)
#' @param times A numeric vector of time points corresponding to columns of `sp_matrix`
#'
#' @return A named vector with "IAE" and "ISE"
#' @export
IAEISE_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (!is.matrix(sp_matrix) && !is.data.frame(sp_matrix)) stop("sp_matrix must be a matrix or data.frame")
  if (length(times) != ncol(sp_matrix)) stop("Length of `times` must match number of columns in `sp_matrix`")
  if (any(diff(times) <= 0)) stop("Time points must be strictly increasing")

  # Average predicted survival curve
  mean_pred_surv <- colMeans(sp_matrix)

  # Empirical KM curve
  km_fit <- survfit(object ~ 1)
  km_data <- data.frame(time = km_fit$time, surv = km_fit$surv)
  km_data <- km_data[km_fit$n.event > 0, , drop = FALSE]

  # Interpolate predicted survival at observed event times
  pred_curve <- data.frame(time = times, surv = mean_pred_surv)
  pred_at_km <- approx(x = pred_curve$time, y = pred_curve$surv,
                       xout = km_data$time, method = "constant", rule = 2)$y

  # Calculate IAE and ISE via trapezoidal approximation
  delta_t <- diff(km_data$time)
  t_IAE <- abs(km_data$surv - pred_at_km)
  t_ISE <- (km_data$surv - pred_at_km)^2

  IAE <- sum(t_IAE[-length(t_IAE)] * delta_t)
  ISE <- sum(t_ISE[-length(t_ISE)] * delta_t)

  result <- c(IAE = round(IAE, 4), ISE = round(ISE, 4))
  return(result)
}



library(survival)
data(sTRACE)

mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)

# Predict survival probabilities at 3 time points
time_points <- c(1, 3, 5)
pred_mat <- predict_aareg(mod, newdata = sTRACE, times = time_points)

# Evaluate IAE and ISE
IAEISE_survmat(Surv(sTRACE$time, sTRACE$status == 9), pred_mat, times = time_points)
