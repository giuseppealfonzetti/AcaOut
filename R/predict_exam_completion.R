#' Predict exam completion probabilities within a horizon
#'
#' For every student-exam pair still pending at the end of `CUTOFF_YEAR`,
#' compute the probability that the exam will be passed by the end of academic
#' year `CUTOFF_YEAR + HORIZON`, conditioning on the truncated history. The
#' forecast mirrors the internal machinery of [predict_k_years()] but keeps the
#' probabilities at the exam level.
#' @inheritParams predict_k_years
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom purrr map_dfr
#' @return A tibble with one row per subject-exam containing the probability of
#'   completing that exam within the target horizon. Exams already completed
#'   before the cutoff are omitted.
#'
#' @export
predict_exam_completion <- function(
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

  cutoff_day <- if (CUTOFF_YEAR == 0L) {
    0L
  } else {
    as.integer(round(CUTOFF_YEAR * YEAR_LENGTH))
  }
  target_year <- CUTOFF_YEAR + HORIZON
  target_day <- as.integer(round(
    (target_year - 1L) * YEAR_LENGTH + GRAD_EXTENSION * YEAR_LENGTH
  ))

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

  yle_trunc <- rep.int(100L, nrow(grades))
  for (i in seq_len(nrow(grades))) {
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
    mu_ability <- rep(0, nrow(covariates))
    mu_speed <- rep(0, nrow(covariates))
  }

  pending_mask <- todo & is.na(grades_trunc)
  subject_id <- FIT$data$labs$obs
  exam_labels <- FIT$data$labs$exams

  out_list <- vector("list", length = nrow(grades))
  for (i in seq_len(nrow(grades))) {
    pending_idx <- which(pending_mask[i, ])
    if (length(pending_idx) == 0L) {
      next
    }
    Ew_row <- Ew[i, ]
    cov_i <- if (n_cov > 0) as.numeric(covariates[i, ]) else numeric(0)

    exam_tbl <- purrr::map_dfr(
      pending_idx,
      function(exam_j) {
        prob_complete <- 0
        for (q in seq_len(nq)) {
          w <- Ew_row[q]
          if (w <= 0 || !is.finite(w)) {
            next
          }
          ability <- nodes[q, 1] + mu_ability[i]
          speed <- nodes[q, 2] + mu_speed[i]

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
          cdf_target <- cpp_pTimeExam(
            EXAM = exam_j - 1L,
            DAY = target_day,
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
          cdf_target <- pmin(pmax(cdf_target, cdf_cut), 1)
          S_cut <- 1 - p_pass * cdf_cut
          S_target <- 1 - p_pass * cdf_target
          if (S_cut <= 1e-12) {
            S_cut <- 1e-12
          }
          if (S_target <= 1e-12) {
            S_target <- 1e-12
          }
          prob_exam <- 1 - (S_target / S_cut)
          prob_exam <- pmin(pmax(prob_exam, 0), 1)
          prob_complete <- prob_complete + w * prob_exam
        }

        tibble::tibble(
          subject_id = subject_id[i],
          cutoff_year = CUTOFF_YEAR,
          horizon = HORIZON,
          target_year = target_year,
          exam = exam_labels[exam_j],
          predicted_completion_prob = prob_complete
        )
      }
    )
    out_list[[i]] <- exam_tbl
  }

  out_tbl <- dplyr::bind_rows(out_list)
  if (nrow(out_tbl) == 0) {
    warning("No pending exams at the selected cutoff; returning empty tibble.")
  }
  out_tbl
}
