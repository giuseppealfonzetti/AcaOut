#' Outcome probabilities using full exam histories
#'
#' Compute the probability that each student experiences graduation, dropout,
#' transfer, or remains enrolled by the end of selected academic years,
#' conditioning on the complete exam history (ignoring observed outcomes). The
#' function mirrors the competing-risk component used during fitting and is
#' useful to assess the calibration of the hazard parameters in isolation.
#'
#' @inheritParams predict_k_years
#' @param TARGET_YEARS Integer vector of academic years for which cumulative
#'   probabilities should be returned. Defaults to the model horizon
#'   (`FIT$data$data_dims$yb`).
#'
#' @return A tibble with one row per student and target year containing
#'   `prob_graduation`, `prob_dropout`, `prob_transfer`, and
#'   `prob_still_enrolled`, together with the posterior EAP ability/speed
#'   computed from exams only.
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom purrr map_dfr
predict_outcome_full <- function(
  FIT,
  TARGET_YEARS = FIT$data$data_dims$yb,
  YEAR_LENGTH = 365,
  GRAD_EXTENSION = 1.5
) {
  if (!(FIT$mod %in% "full")) {
    stop("Outcome prediction is only available for `mod = \"full\"` fits.")
  }
  target_years <- sort(unique(as.integer(TARGET_YEARS)))
  if (any(!is.finite(target_years)) || length(target_years) == 0) {
    stop("`TARGET_YEARS` must contain finite integers.")
  }
  if (any(target_years <= 0)) {
    stop("`TARGET_YEARS` must be positive integers.")
  }
  if (length(YEAR_LENGTH) != 1 || !is.finite(YEAR_LENGTH) || YEAR_LENGTH <= 0) {
    stop("`YEAR_LENGTH` must be a positive scalar.")
  }
  if (length(GRAD_EXTENSION) != 1 || !is.finite(GRAD_EXTENSION) || GRAD_EXTENSION <= 1) {
    stop("`GRAD_EXTENSION` must be greater than 1.")
  }

  data_dims <- FIT$data$data_dims
  if (any(target_years > data_dims$yb)) {
    stop("`TARGET_YEARS` cannot exceed the model horizon (", data_dims$yb, ").")
  }

  grades <- FIT$data$gradesMat
  times <- FIT$data$timeMat
  todo <- FIT$data$todoMat
  covariates <- as.matrix(FIT$data$X)
  max_time_full <- FIT$data$max_time
  n_students <- nrow(grades)

  # Treat outcomes as unobserved so that predictions depend only on exams
  outcome_null <- rep.int(0L, n_students)
  last_year_null <- rep.int(0L, n_students)

  yle <- rep.int(100L, n_students)
  for (i in seq_len(n_students)) {
    required <- todo[i, ]
    obs_times <- times[i, required]
    obs_times <- obs_times[is.finite(obs_times)]
    if (length(obs_times) > 0) {
      yle[i] <- as.integer(ceiling(max(obs_times) / YEAR_LENGTH))
    }
  }

  max_time_used <- max_time_full
  max_time_used[!is.finite(max_time_used)] <- 0L

  estep_out <- cpp_estep(
    THETA = FIT$fit$par,
    EXAMS_GRADES = grades,
    EXAMS_DAYS = times,
    EXAMS_SET = todo,
    EXAMS_OBSFLAG = !is.na(times),
    MAX_DAY = max_time_used,
    OUTCOME = outcome_null,
    EXT_COVARIATES = covariates,
    YEAR_FIRST = FIT$data$first_year,
    YEAR_LAST = last_year_null,
    YEAR_LAST_EXAM = yle,
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

  target_days <- as.integer(round((target_years - 1L) * YEAR_LENGTH +
    GRAD_EXTENSION * YEAR_LENGTH))
  subject_id <- FIT$data$labs$obs

  purrr::map_dfr(seq_len(n_students), function(i) {
    Ew_row <- Ew[i, ]
    cov_i <- if (n_cov > 0) as.numeric(covariates[i, ]) else numeric(0)
    pending_idx <- which(todo[i, ] & is.na(grades[i, ]))
    finished <- length(pending_idx) == 0L
    cutoff_day <- max_time_used[i]

    ability_eap <- sum(Ew_row * (nodes[, 1] + mu_ability[i]))
    speed_eap <- sum(Ew_row * (nodes[, 2] + mu_speed[i]))

    prob_grad_years <- prob_drop_years <- prob_trans_years <- numeric(length(target_years))

    for (q in seq_len(nq)) {
      w <- Ew_row[q]
      if (w <= 0 || !is.finite(w)) next
      ability <- nodes[q, 1] + mu_ability[i]
      speed <- nodes[q, 2] + mu_speed[i]

      hazard_grad_vec <- hazard_drop_vec <- hazard_trans_vec <- numeric(length(target_years))
      for (t_idx in seq_along(target_years)) {
        target_year <- target_years[t_idx]
        hazard_grad_vec[t_idx] <- cpp_hazard(
          1L, target_year, theta, cov_i, ability, speed, data_dims$yb, FALSE
        )$prob
        hazard_drop_vec[t_idx] <- cpp_hazard(
          2L, target_year, theta, cov_i, ability, speed, data_dims$yb, FALSE
        )$prob
        hazard_trans_vec[t_idx] <- cpp_hazard(
          3L, target_year, theta, cov_i, ability, speed, data_dims$yb, FALSE
        )$prob
      }

      if (finished) {
        prob_grad_years <- prob_grad_years + w * hazard_grad_vec
        next
      }

      prob_finish_vec <- rep(1, length(target_years))
      for (exam_j in pending_idx) {
        p_pass <- cpp_pGreaterGrades(
          1L, exam_j - 1L, theta, cov_i, data_dims$n_grades, data_dims$n_exams, ability, FALSE, FALSE
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
        if (S_cut <= 1e-12) S_cut <- 1e-12
        for (t_idx in seq_along(target_days)) {
          cdf_target <- cpp_pTimeExam(
            EXAM = exam_j - 1L,
            DAY = target_days[t_idx],
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
          if (S_target <= 1e-12) S_target <- 1e-12
          prob_exam <- 1 - (S_target / S_cut)
          prob_exam <- min(max(prob_exam, 0), 1)
          prob_finish_vec[t_idx] <- prob_finish_vec[t_idx] * prob_exam
        }
      }

      prob_finish_vec <- pmin(pmax(prob_finish_vec, 0), 1)
      prob_grad_years <- prob_grad_years + w * (prob_finish_vec * hazard_grad_vec)
      drop_contrib <- (1 - prob_finish_vec) * hazard_drop_vec
      trans_contrib <- (1 - prob_finish_vec) * hazard_trans_vec
      prob_drop_years <- prob_drop_years + w * drop_contrib
      prob_trans_years <- prob_trans_years + w * trans_contrib
    }

    prob_grad_years <- pmin(pmax(prob_grad_years, 0), 1)
    prob_drop_years <- pmin(pmax(prob_drop_years, 0), 1)
    prob_trans_years <- pmin(pmax(prob_trans_years, 0), 1)

    survival <- 1
    cum_grad <- cum_drop <- cum_trans <- 0
    prob_grad_cum <- prob_drop_cum <- prob_trans_cum <- prob_still <- numeric(length(target_years))

    for (t_idx in seq_along(target_years)) {
      g <- prob_grad_years[t_idx]
      d <- prob_drop_years[t_idx]
      tr <- prob_trans_years[t_idx]
      total <- pmin(pmax(g + d + tr, 0), 1)

      cum_grad <- cum_grad + survival * g
      cum_drop <- cum_drop + survival * d
      cum_trans <- cum_trans + survival * tr
      survival <- survival * (1 - total)
      survival <- pmin(pmax(survival, 0), 1)

      prob_grad_cum[t_idx] <- cum_grad
      prob_drop_cum[t_idx] <- cum_drop
      prob_trans_cum[t_idx] <- cum_trans
      prob_still[t_idx] <- survival
    }

    tibble(
      subject_id = subject_id[i],
      target_year = target_years,
      prob_graduation = prob_grad_cum,
      prob_dropout = prob_drop_cum,
      prob_transfer = prob_trans_cum,
      prob_still_enrolled = prob_still,
      ability_eap = ability_eap,
      speed_eap = speed_eap
    )
  })
}
