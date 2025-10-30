#' K-year-ahead outcome predictions
#'
#' Compute the probability that each student experiences graduation, dropout,
#' transfer, or remains enrolled by the end of academic year
#' `CUTOFF_YEAR + HORIZON`, conditioning on the information available up to the
#' end of `CUTOFF_YEAR`. The function reuses the posterior weights computed on
#' the truncated history (as in [predict_one_year()]) and propagates the hazards
#' across multiple future years without leaking information from later periods.
#'
#' @inheritParams predict_one_year
#' @param HORIZON Strictly positive integer specifying how many academic years
#'   ahead to evaluate. The final horizon year equals `CUTOFF_YEAR + HORIZON`.
#'
#' @return A tibble with one row per student containing cumulative probabilities
#'   of graduation, dropout, transfer, and remaining enrolled by the end of the
#'   target horizon.
#'
#' @export
#'
#' @importFrom tibble tibble
predict_k_years <- function(
  FIT,
  CUTOFF_YEAR,
  HORIZON,
  YEAR_LENGTH = 365,
  GRAD_EXTENSION = 1.5
) {
  if (!(FIT$mod %in% "full")) {
    stop("Prediction is only available for `mod = \"full\"` fits.")
  }
  if (length(CUTOFF_YEAR) != 1 || !is.finite(CUTOFF_YEAR)) {
    stop("`CUTOFF_YEAR` must be a single finite value.")
  }
  CUTOFF_YEAR <- as.integer(CUTOFF_YEAR)
  if (CUTOFF_YEAR < 0) {
    stop("`CUTOFF_YEAR` must be non-negative.")
  }
  if (length(HORIZON) != 1 || !is.finite(HORIZON)) {
    stop("`HORIZON` must be a single finite value.")
  }
  HORIZON <- as.integer(HORIZON)
  if (HORIZON <= 0) {
    stop("`HORIZON` must be strictly positive.")
  }
  if (length(YEAR_LENGTH) != 1 || !is.finite(YEAR_LENGTH) || YEAR_LENGTH <= 0) {
    stop("`YEAR_LENGTH` must be a positive scalar.")
  }
  if (
    length(GRAD_EXTENSION) != 1 ||
      !is.finite(GRAD_EXTENSION) ||
      GRAD_EXTENSION <= 1
  ) {
    stop("`GRAD_EXTENSION` must be greater than 1.")
  }

  data_dims <- FIT$data$data_dims
  yb <- data_dims$yb
  if (CUTOFF_YEAR + HORIZON > yb) {
    stop("`CUTOFF_YEAR + HORIZON` exceeds the maximum modelled year.")
  }

  subject_id <- FIT$data$labs$obs
  n_students <- length(subject_id)
  target_years <- CUTOFF_YEAR + seq_len(HORIZON)
  target_days <- as.integer(round(
    (target_years - 1L) * YEAR_LENGTH + GRAD_EXTENSION * YEAR_LENGTH
  ))
  cutoff_day <- if (CUTOFF_YEAR == 0L) {
    0L
  } else {
    as.integer(round(CUTOFF_YEAR * YEAR_LENGTH))
  }

  grades <- FIT$data$gradesMat
  times <- FIT$data$timeMat
  todo <- FIT$data$todoMat
  outcome_full <- FIT$data$outcome
  last_year_full <- FIT$data$last_year
  max_time_full <- FIT$data$max_time
  covariates <- as.matrix(FIT$data$X)

  grades_trunc <- grades
  times_trunc <- times
  mask_future <- !is.na(times_trunc) & times_trunc > cutoff_day
  grades_trunc[mask_future] <- NA_integer_
  times_trunc[mask_future] <- NA_integer_

  exam_limit <- if (cutoff_day > 0L) cutoff_day else 1L
  max_time_trunc <- as.integer(pmax(pmin(max_time_full, exam_limit), 1L))
  last_year_trunc <- as.integer(pmin(last_year_full, CUTOFF_YEAR))
  outcome_trunc <- ifelse(last_year_full > CUTOFF_YEAR, 0L, outcome_full)

  yle_trunc <- rep.int(100L, n_students)
  for (i in seq_len(n_students)) {
    required <- todo[i, ]
    obs_times <- times_trunc[i, required]
    if (length(obs_times) > 0 && all(!is.na(obs_times))) {
      yle_trunc[i] <- as.integer(ceiling(max(obs_times) / YEAR_LENGTH))
    }
  }

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

  theta <- FIT$fit$par
  dim_irt <- FIT$data$par_dims$grtcm
  n_cov <- data_dims$n_cov
  if (n_cov > 0) {
    gamma_ability <- theta[(dim_irt + 3):(dim_irt + 2 + n_cov)]
    gamma_speed <- theta[(dim_irt + 3 + n_cov):(dim_irt + 2 + 2 * n_cov)]
    mu_ability <- as.numeric(covariates %*% gamma_ability)
    mu_speed <- as.numeric(covariates %*% gamma_speed)
  } else {
    mu_ability <- rep(0, n_students)
    mu_speed <- rep(0, n_students)
  }

  pending <- todo & is.na(grades_trunc)
  horizon_len <- length(target_years)
  prob_grad_years <- matrix(0, n_students, horizon_len)
  prob_drop_years <- matrix(0, n_students, horizon_len)
  prob_trans_years <- matrix(0, n_students, horizon_len)

  for (i in seq_len(n_students)) {
    cov_i <- if (n_cov > 0) as.numeric(covariates[i, ]) else numeric(0)
    pending_idx <- which(pending[i, ])
    finished <- length(pending_idx) == 0L
    Ew_row <- Ew[i, ]

    for (q in seq_len(nq)) {
      w <- Ew_row[q]
      if (w <= 0 || !is.finite(w)) {
        next
      }
      ability <- nodes[q, 1] + mu_ability[i]
      speed <- nodes[q, 2] + mu_speed[i]

      hazard_grad_vec <- numeric(horizon_len)
      hazard_drop_vec <- numeric(horizon_len)
      hazard_trans_vec <- numeric(horizon_len)
      for (t in seq_len(horizon_len)) {
        target_year <- target_years[t]
        hazard_grad_vec[t] <- cpp_hazard(
          1L,
          target_year,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
        hazard_drop_vec[t] <- cpp_hazard(
          2L,
          target_year,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
        hazard_trans_vec[t] <- cpp_hazard(
          3L,
          target_year,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
      }

      if (finished) {
        prob_grad_years[i, ] <- prob_grad_years[i, ] + w * hazard_grad_vec
        next
      }

      prob_finish_vec <- rep(1, horizon_len)
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
        cdf_cut <- cpp_pTimeExam(
          EXAM = exam_j - 1L,
          DAY = cutoff_day,
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
        cdf_cut <- pmin(pmax(cdf_cut, 0), 1)
        S_cut <- 1 - p_pass * cdf_cut
        if (S_cut <= 1e-12) {
          S_cut <- 1e-12
        }
        for (t in seq_len(horizon_len)) {
          cdf_target <- cpp_pTimeExam(
            EXAM = exam_j - 1L,
            DAY = target_days[t],
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
          cdf_target <- pmin(pmax(cdf_target, cdf_cut), 1)
          S_target <- 1 - p_pass * cdf_target
          if (S_target <= 1e-12) {
            S_target <- 1e-12
          }
          prob_exam <- 1 - (S_target / S_cut)
          prob_exam <- min(max(prob_exam, 0), 1)
          prob_finish_vec[t] <- prob_finish_vec[t] * prob_exam
        }
      }
      prob_finish_vec <- pmin(pmax(prob_finish_vec, 0), 1)

      prob_grad_years[i, ] <- prob_grad_years[i, ] +
        w * (prob_finish_vec * hazard_grad_vec)
      drop_contrib <- (1 - prob_finish_vec) * hazard_drop_vec
      trans_contrib <- (1 - prob_finish_vec) * hazard_trans_vec
      prob_drop_years[i, ] <- prob_drop_years[i, ] + w * drop_contrib
      prob_trans_years[i, ] <- prob_trans_years[i, ] + w * trans_contrib
    }
  }

  prob_grad_years <- pmin(pmax(prob_grad_years, 0), 1)
  prob_drop_years <- pmin(pmax(prob_drop_years, 0), 1)
  prob_trans_years <- pmin(pmax(prob_trans_years, 0), 1)

  survival <- rep(1, n_students)
  cum_grad <- cum_drop <- cum_trans <- rep(0, n_students)
  for (t in seq_len(horizon_len)) {
    g <- prob_grad_years[, t]
    d <- prob_drop_years[, t]
    tr <- prob_trans_years[, t]
    total <- pmin(pmax(g + d + tr, 0), 1)

    cum_grad <- cum_grad + survival * g
    cum_drop <- cum_drop + survival * d
    cum_trans <- cum_trans + survival * tr
    survival <- survival * (1 - total)
    survival <- pmin(pmax(survival, 0), 1)
  }

  prob_still <- survival

  occurred_before_cutoff <- outcome_full %in%
    c(1L, 2L, 3L) &
    last_year_full <= CUTOFF_YEAR
  if (any(occurred_before_cutoff)) {
    grad_idx <- occurred_before_cutoff & outcome_full == 1L
    drop_idx <- occurred_before_cutoff & outcome_full == 2L
    trans_idx <- occurred_before_cutoff & outcome_full == 3L
    cum_grad[grad_idx] <- 1
    cum_drop[drop_idx] <- 1
    cum_trans[trans_idx] <- 1
    prob_still[occurred_before_cutoff] <- 0
  }

  tibble(
    subject_id = subject_id,
    cutoff_year = CUTOFF_YEAR,
    horizon = HORIZON,
    target_year = CUTOFF_YEAR + HORIZON,
    n_exams_done = rowSums(!is.na(times_trunc)),
    n_exams_todo = rowSums(todo),
    prob_graduation = pmin(pmax(cum_grad, 0), 1),
    prob_dropout = pmin(pmax(cum_drop, 0), 1),
    prob_transfer = pmin(pmax(cum_trans, 0), 1),
    prob_still_enrolled = pmin(pmax(prob_still, 0), 1)
  )
}
