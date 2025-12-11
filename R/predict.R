#' Predicted outcome probabilities one year ahead
#'
#' @param FIT Fitted object returned by [fit_EM()] or [fit_BFGS()].
#' @param CUTOFF_YEAR Consider data up to the end of `CUTOFF_YEAR`.
#' @param YEAR_LENGTH Number of day within each year. Censoring day is computed as `CUTOFF_YEAR*YEAR_LENGTH`.
#' @param GRAD_EXTENSION Scalar multiplier to account for extra graduation sessions. Last day useful for graduation is computed as `CUTOFF_YEAR*YEAR_LENGTH+ GRAD_EXTENSION*YEAR_LENGTH`.
#'
#' @export
predict_one_year <- function(
  FIT,
  CUTOFF_YEAR,
  YEAR_LENGTH = 365,
  GRAD_EXTENSION = 1.5
) {
  # input checks
  if (!(FIT$mod %in% "full")) {
    stop("Prediction is only available for `mod = \"full\"` fits.")
  }
  if (length(CUTOFF_YEAR) != 1 || !is.finite(CUTOFF_YEAR)) {
    stop(
      "`CUTOFF_YEAR` must be a single finite value."
    )
  }
  CUTOFF_YEAR <- as.integer(CUTOFF_YEAR)
  if (CUTOFF_YEAR < 0) {
    stop("`CUTOFF_YEAR` must be non-negative.")
  }
  if (length(YEAR_LENGTH) != 1 || !is.finite(YEAR_LENGTH) || YEAR_LENGTH <= 0) {
    stop(
      "`YEAR_LENGTH` must be a positive scalar."
    )
  }
  if (
    length(GRAD_EXTENSION) != 1 ||
      !is.finite(GRAD_EXTENSION) ||
      GRAD_EXTENSION <= 1
  ) {
    stop("`GRAD_EXTENSION` must be greater than 1.")
  }

  # read data
  grades <- FIT$data$gradesMat
  times <- FIT$data$timeMat
  todo <- FIT$data$todoMat
  outcome_full <- FIT$data$outcome
  last_year_full <- FIT$data$last_year
  max_time_full <- FIT$data$max_time
  covariates <- as.matrix(FIT$data$X)
  labs <- FIT$data$labs
  data_dims <- FIT$data$data_dims
  par_dims <- FIT$data$par_dims

  # set up cutoff and truncate
  cutoff_exam_day <- if (CUTOFF_YEAR == 0L) {
    0L
  } else {
    as.integer(round(
      GRAD_EXTENSION * YEAR_LENGTH + (CUTOFF_YEAR - 1) * YEAR_LENGTH
    ))
  }
  grades_trunc <- grades
  times_trunc <- times
  mask_future <- !is.na(times_trunc) & times_trunc > cutoff_exam_day
  grades_trunc[mask_future] <- NA
  times_trunc[mask_future] <- NA
  exam_limit <- if (cutoff_exam_day > 0L) cutoff_exam_day else 1L
  max_time_trunc <- pmin(max_time_full, exam_limit)
  max_time_trunc <- pmax(max_time_trunc, 1L)
  last_year_trunc <- pmin(last_year_full, CUTOFF_YEAR)
  outcome_trunc <- ifelse(last_year_full > CUTOFF_YEAR, 0L, outcome_full)
  yle_trunc <- rep(100L, nrow(grades))
  for (i in seq_len(nrow(grades))) {
    required <- todo[i, ]
    obs_times <- times_trunc[i, required]
    if (all(!is.na(obs_times)) && length(obs_times) > 0) {
      yle_trunc[i] <- as.integer(ceiling(max(obs_times) / YEAR_LENGTH))
    }
  }

  # Compute integration weights
  estep_out <- cpp_estep(
    THETA = FIT$fit$par,
    EXAMS_GRADES = grades_trunc,
    EXAMS_DAYS = times_trunc,
    EXAMS_SET = todo,
    EXAMS_OBSFLAG = !is.na(times_trunc),
    MAX_DAY = max_time_trunc,
    OUTCOME = outcome_trunc,
    EXT_COVARIATES = covariates,
    YEAR_FIRST = FIT$data$first_year,
    YEAR_LAST = last_year_trunc,
    YEAR_LAST_EXAM = yle_trunc,
    GRID = FIT$grid,
    WEIGHTS = FIT$weights,
    YB = data_dims$yb,
    N_GRADES = data_dims$n_grades,
    N_EXAMS = data_dims$n_exams,
    MOD = FIT$mod
  )
  Ew <- estep_out$Ew
  nodes <- estep_out$nodes
  nq <- nrow(nodes)
  Ew[!is.finite(Ew)] <- 0
  row_sums <- rowSums(Ew)
  zero_rows <- row_sums <= 0
  if (any(zero_rows)) {
    Ew[zero_rows, ] <- 1 / nq
    row_sums[zero_rows] <- 1
  }
  Ew <- Ew / row_sums

  # set up for computing outcome probs for each student
  theta <- FIT$fit$par
  dim_irt <- par_dims$grtcm
  n_cov <- data_dims$n_cov
  if (n_cov > 0) {
    gamma_ability <- theta[(dim_irt + 3):(dim_irt + 2 + n_cov)]
    gamma_speed <- theta[(dim_irt + 3 + n_cov):(dim_irt + 2 + 2 * n_cov)]
    mu_ability <- as.numeric(covariates %*% gamma_ability)
    mu_speed <- as.numeric(covariates %*% gamma_speed)
  } else {
    mu_ability <- rep(0, nrow(covariates))
    mu_speed <- rep(0, nrow(covariates))
  }
  target_year <- CUTOFF_YEAR + 1L
  grad_end_day <- as.integer(round(
    CUTOFF_YEAR * YEAR_LENGTH + GRAD_EXTENSION * YEAR_LENGTH
  ))
  prob_mat <- matrix(
    0,
    nrow = nrow(covariates),
    ncol = 3,
    dimnames = list(NULL, c("graduation", "dropout", "transfer"))
  )
  pending <- todo & is.na(grades_trunc)

  # compute probability mat
  for (i in seq_len(nrow(covariates))) {
    cov_i <- if (n_cov > 0) as.numeric(covariates[i, ]) else numeric(0)
    pending_idx <- which(pending[i, ])
    finished <- length(pending_idx) == 0L

    # loop over quadrature nodes
    for (q in seq_len(nq)) {
      w <- Ew[i, q]
      if (w <= 0 || !is.finite(w)) {
        next
      }
      ability <- nodes[q, 1] + mu_ability[i]
      speed <- nodes[q, 2] + mu_speed[i]
      hazard_g <- cpp_hazard(
        1L,
        target_year,
        theta,
        cov_i,
        ability,
        speed,
        data_dims$yb,
        FALSE
      )$prob
      if (finished) {
        prob_grad <- hazard_g
        prob_drop <- 0
        prob_trans <- 0
      } else {
        hazard_d <- cpp_hazard(
          2L,
          target_year,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
        hazard_t <- cpp_hazard(
          3L,
          target_year,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
        prob_finish <- 1
        for (exam_j in pending_idx) {
          p_pass <- cpp_pGreaterGrades(
            1L,
            exam_j - 1L,
            theta,
            cov_i,
            data_dims$n_grades,
            data_dims$n_exams,
            ability,
            FALSE,
            FALSE
          )$prob
          cdf_y <- cpp_pTimeExam(
            EXAM = exam_j - 1L,
            DAY = cutoff_exam_day,
            THETA = theta,
            COVARIATES = cov_i,
            N_GRADES = data_dims$n_grades,
            N_EXAMS = data_dims$n_exams,
            SPEED = speed,
            ABILITY = ability,
            CDFFLAG = TRUE,
            LOGFLAG = FALSE,
            LATPARFLAG = FALSE
          )$prob
          cdf_yplus <- cpp_pTimeExam(
            EXAM = exam_j - 1L,
            DAY = grad_end_day,
            THETA = theta,
            COVARIATES = cov_i,
            N_GRADES = data_dims$n_grades,
            N_EXAMS = data_dims$n_exams,
            SPEED = speed,
            ABILITY = ability,
            CDFFLAG = TRUE,
            LOGFLAG = FALSE,
            LATPARFLAG = FALSE
          )$prob
          cdf_y <- pmin(pmax(cdf_y, 0), 1)
          cdf_yplus <- pmin(pmax(cdf_yplus, 0), 1)
          S_y <- 1 - p_pass * cdf_y
          if (S_y <= 1e-12) {
            S_y <- 1e-12
          }
          prob_exam <- 1 - (1 - p_pass * cdf_yplus) / S_y
          prob_exam <- pmin(pmax(prob_exam, 0), 1)
          prob_finish <- prob_finish * prob_exam
        }
        prob_grad <- prob_finish * hazard_g
        prob_drop <- (1 - prob_finish) * hazard_d
        prob_trans <- (1 - prob_finish) * hazard_t
      }
      prob_mat[i, ] <- prob_mat[i, ] + w * c(prob_grad, prob_drop, prob_trans)
    }
  }

  # correct for numerical stability
  prob_mat <- pmin(pmax(prob_mat, 0), 1)
  prob_still <- pmax(0, 1 - rowSums(prob_mat))

  # check for outcome already occurrd
  occurred <- last_year_full <= CUTOFF_YEAR & outcome_full %in% c(1L, 2L, 3L)
  for (idx in which(occurred)) {
    prob_mat[idx, ] <- 0
    prob_mat[idx, outcome_full[idx]] <- 0
    prob_still[idx] <- 0
  }

  # return the mat
  out <- tibble::tibble(
    subject_id = labs$obs,
    CUTOFF_YEAR = CUTOFF_YEAR,
    target_year = target_year,
    outcome = factor(
      outcome_full,
      levels = 0:3,
      labels = c("enrolled", "graduated", "dropout", "transfer")
    ),
    last_year = last_year_full,
    prob_graduation = prob_mat[, "graduation"],
    prob_dropout = prob_mat[, "dropout"],
    prob_transfer = prob_mat[, "transfer"],
    prob_still_enrolled = prob_still,
    event_in_target_year = dplyr::case_when(
      outcome == "graduated" & last_year == target_year ~ "graduated",
      outcome == "dropout" & last_year == target_year ~ "dropout",
      outcome == "transfer" & last_year == target_year ~ "transfer",
      outcome == "graduated" & last_year < target_year ~ "already_graduated",
      outcome == "dropout" & last_year < target_year ~ "already_dropout",
      outcome == "transfer" & last_year < target_year ~ "already_transfer",
      TRUE ~ "not_yet"
    )
  )

  return(out)
}
