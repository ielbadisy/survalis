library(data.table)

Cindex_survmat <- function(object, predicted, t_star = NULL) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")

  time <- object[, 1]
  status <- object[, 2]

  surv_prob <- if (!is.null(t_star)) {
    t_name <- paste0("t=", t_star)
    if (!(t_name %in% colnames(predicted))) stop("t_star = ", t_star, " not found.")
    predicted[[t_name]]
  } else {
    predicted[[ncol(predicted)]]
  }

  risk_score <- 1 - surv_prob
  dt <- data.table(id = seq_along(time), time, status, risk_score)

  #self-join for pairwise comparisons
  dt_pairs <- dt[dt, on = .(id < id), allow.cartesian = TRUE]

  dt_pairs <- dt_pairs[
    !((time < i.time & status == 0) | (i.time < time & i.status == 0)) &  # both uncensored or valid
      !(time == i.time & status == 0 & i.status == 0)
  ]

  dt_pairs[, comp := fifelse(
    time != i.time,
    fifelse((time < i.time & risk_score > i.risk_score) |
              (i.time < time & i.risk_score > risk_score), 1,
      fifelse(risk_score == i.risk_score, 0.5, 0)),
    fifelse(status + i.status > 0,
      fifelse(risk_score == i.risk_score, 1, 0.5),
      NA_real_)
  )]

  permissible <- nrow(dt_pairs[!is.na(comp)])
  concord <- sum(dt_pairs[!is.na(comp), comp])
  C_index <- concord / permissible

  names(C_index) <- "C index"
  return(round(C_index, 6))
}

mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)
pred <- predict_aareg(mod, newdata = sTRACE, times = c(1, 3, 5))


Brier <- function(object, pre_sp, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(pre_sp) != nrow(object)) stop("Length of predictions must match number of observations")

  dt <- data.table(time = object[, 1], status = object[, 2], sp = pre_sp)

  # estimate censoring distribution G(t) using Kaplan-Meier on 1 - status
  km_fit <- survfit(Surv(time, 1 - status) ~ 1, data = dt)
  G_vals <- summary(km_fit, times = sort(unique(c(t_star, dt$time))), extend = TRUE)$surv
  names(G_vals) <- as.character(sort(unique(c(t_star, dt$time))))

  Gt <- G_vals[as.character(t_star)]
  if (is.na(Gt) || Gt == 0) stop("Gt(t_star) is NA or zero")

  dt[, Gti := G_vals[as.character(time)]]
  dt[, brier := 0]

  # event before t_star: use pre_sp^2 / G(ti)
  dt[time < t_star & status == 1 & !is.na(Gti) & Gti > 0,
     brier := (sp^2) / Gti]

  # still at risk at t_star
  dt[time >= t_star & Gt > 0, brier := ((1 - sp)^2) / Gt]

  bs <- mean(dt$brier, na.rm = TRUE)
  names(bs) <- "Brier Score"
  round(bs, 6)
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
  IBS <- integral / (max(times) - min(times))
  names(IBS) <- "IBS"
  round(IBS, 6)
}

IAEISE_survmat <- function(object, sp_matrix, times) {
  if (!inherits(object, "Surv")) stop("object must be a survival object")
  if (length(times) != ncol(sp_matrix)) stop("Length of times must match sp_matrix columns")

  mean_pred_surv <- colMeans(sp_matrix)
  km_fit <- survfit(object ~ 1)
  km_dt <- data.table(time = km_fit$time, surv = km_fit$surv)[km_fit$n.event > 0]

  # interpolate predicted survival at km times
  km_dt[, pred := approx(x = times, y = mean_pred_surv, xout = time, method = "constant", rule = 2)$y]
  km_dt[, `:=`(
    iae = abs(surv - pred),
    ise = (surv - pred)^2
  )]

  # integrate over time
  km_dt[, dt := shift(time, type = "lead") - time]
  iae_val <- km_dt[!is.na(dt), sum(iae * dt)]
  ise_val <- km_dt[!is.na(dt), sum(ise * dt)]

  round(c(IAE = iae_val, ISE = ise_val), 4)
}


mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf, data = sTRACE)
times <- c(1, 2, 5)
evaluate_survlearner(mod, metrics = c("cindex", "ibs", "iae", "ise"), times = times)
