#' Concordance index for right-censored survival based on survival probability matrix
#'
#' @param object A survival object created with `Surv(time, status)`
#' @param predicted A data.frame or matrix of survival probabilities (output from `predict_aareg`)
#' @param t_star Time point at which to evaluate the C-index (must match column name format `t=3`)
#'
#' @return Estimated concordance index (C-index)
#' @export
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
      surv_prob <- predicted[[ncol(predicted)]]  # fallback: last time column
      warning("No t_star provided, using last column: ", colnames(predicted)[ncol(predicted)])
  }
  
  if (length(surv_prob) != length(time)) stop("Length of predictions does not match time vector.")
  
  # Convert survival probability to risk score (higher = more risk)
  risk_score <- 1 - surv_prob
  
  # --- Concordance logic ---
  permissible <- 0
  concord <- 0
  par_concord <- 0
  
  n <- length(time)
  for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
          if ((time[i] < time[j] & status[i] == 0) | (time[j] < time[i] & status[j] == 0)) {
              next
          }
          
          if (time[i] == time[j] & status[i] == 0 & status[j] == 0) {
              next
          }
          
          permissible <- permissible + 1
          
          if (time[i] != time[j]) {
              if ((time[i] < time[j] & risk_score[i] > risk_score[j]) |
                  (time[j] < time[i] & risk_score[j] > risk_score[i])) {
                  concord <- concord + 1
              } else if (risk_score[i] == risk_score[j]) {
                  par_concord <- par_concord + 0.5
              }
          }
          
          if (time[i] == time[j] & status[i] == 1 & status[j] == 1) {
              if (risk_score[i] == risk_score[j]) {
                  concord <- concord + 1
              } else {
                  par_concord <- par_concord + 0.5
              }
          }
          
          if (time[i] == time[j] &
              ((status[i] == 1 & status[j] == 0) | (status[i] == 0 & status[j] == 1))) {
              if ((status[i] == 1 & risk_score[i] > risk_score[j]) |
                  (status[j] == 1 & risk_score[j] > risk_score[i])) {
                  concord <- concord + 1
              } else {
                  par_concord <- par_concord + 0.5
              }
          }
      }
  }
  
  C_index <- (concord + par_concord) / permissible
  names(C_index) <- "C index"
  return(round(C_index, 6))
}



library(survival)
library(timereg)

data(sTRACE)

mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf,
                 data = sTRACE, max.time = 7)

pred <- predict_aareg(mod, newdata = sTRACE, times = c(1, 3, 5))

# Compute C-index at t = 3
Cindex_survmat(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3)


#************

# use data.table

#' Fast concordance index for survival probability matrix using data.table
#'
#' @param object A survival object created with `Surv(time, status)`
#' @param predicted A data.frame or matrix of survival probabilities (output from `predict_aareg`)
#' @param t_star Time point at which to evaluate the C-index
#'
#' @return Estimated concordance index (C-index)
#' @import data.table
#' @export
Cindex_survmat_fast <- function(object, predicted, t_star = NULL) {
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
    warning("No t_star provided, using last column: ", colnames(predicted)[ncol(predicted)])
  }
  
  risk <- 1 - surv_prob
  
  dt <- data.table(id = 1:length(risk), time, status, risk)
  
  setkey(dt, time)

  # Cross join all pairs (i, j) with i < j
  dt_pairs <- CJ(i = dt$id, j = dt$id)[i < j]
  dt_pairs <- dt_pairs[
    , `:=`(
      time_i = dt$time[i],
      time_j = dt$time[j],
      status_i = dt$status[i],
      status_j = dt$status[j],
      risk_i = dt$risk[i],
      risk_j = dt$risk[j]
    )
  ]

  # Define permissible pairs
  dt_pairs <- dt_pairs[!((time_i < time_j & status_i == 0) |
                         (time_j < time_i & status_j == 0) |
                         (time_i == time_j & status_i == 0 & status_j == 0))]

  dt_pairs[, comp := 0]

  # Concordant: risk_i > risk_j when i dies before j (or vice versa)
  dt_pairs[time_i < time_j & risk_i > risk_j, comp := 1]
  dt_pairs[time_j < time_i & risk_j > risk_i, comp := 1]

  # Tied risks
  dt_pairs[time_i != time_j & risk_i == risk_j, comp := 0.5]

  # Equal time: both died
  dt_pairs[time_i == time_j & status_i == 1 & status_j == 1,
           comp := ifelse(risk_i == risk_j, 1, 0.5)]

  # Equal time: one died, one censored
  dt_pairs[time_i == time_j &
             ((status_i == 1 & status_j == 0) | (status_i == 0 & status_j == 1)),
           comp := ifelse((status_i == 1 & risk_i > risk_j) |
                          (status_j == 1 & risk_j > risk_i), 1, 0.5)]

  # Output concordance
  concord <- sum(dt_pairs$comp)
  permissible <- nrow(dt_pairs)

  round(concord / permissible, 6)
}


### with lapply


Cindex_apply <- function(object, predicted, t_star = NULL) {
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
    warning("No t_star provided, using last column: ", colnames(predicted)[ncol(predicted)])
  }
  
  risk <- 1 - surv_prob
  n <- length(time)
  
  results <- lapply(1:(n - 1), function(i) {
    sapply((i + 1):n, function(j) {
      # Skip non-permissible pairs
      if ((time[i] < time[j] && status[i] == 0) ||
          (time[j] < time[i] && status[j] == 0) ||
          (time[i] == time[j] && status[i] == 0 && status[j] == 0)) {
        return(NA_real_)
      }

      # Concordance logic
      if (time[i] != time[j]) {
        if ((time[i] < time[j] && risk[i] > risk[j]) ||
            (time[j] < time[i] && risk[j] > risk[i])) {
          return(1)
        } else if (risk[i] == risk[j]) {
          return(0.5)
        } else {
          return(0)
        }
      }

      # Equal times
      if (status[i] + status[j] > 0) {
        if (risk[i] == risk[j]) return(1)
        else return(0.5)
      }

      return(0)
    })
  })

  res_vec <- unlist(results)
  res_vec <- res_vec[!is.na(res_vec)]
  cindex <- mean(res_vec)
  round(cindex, 6)
}



Cindex_survmat(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3)


Cindex_apply(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3)


Cindex_survmat_fast(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3)




library(survival)
library(timereg)
library(data.table)
library(microbenchmark)

# Run benchmark
benchmark <- microbenchmark(
  Cindex_survmat(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3),
  Cindex_apply(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3),
  Cindex_survmat_fast(Surv(sTRACE$time, sTRACE$status == 9), predicted = pred, t_star = 3),
  times = 5L
)

print(benchmark)

