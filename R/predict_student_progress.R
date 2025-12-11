#' Sequential risk updates for a single student
#'
#' @inheritParams predict_k_years
#' @param SUBJECT_ID Character identifier of the student (must match
#'   `FIT$data$labs$obs`).
#' @param horizons Integer vector of academic years to forecast relative to the
#'   current cutoff. By default the function projects to every remaining year.
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows mutate select filter arrange
#' @importFrom purrr map_dbl map map_dfr
predict_student_progress <- function(
  FIT,
  SUBJECT_ID,
  horizons = NULL,
  YEAR_LENGTH = 365,
  GRAD_EXTENSION = 1.5
) {
  if (!(FIT$mod %in% "full")) {
    stop("Sequential prediction is only available for `mod = \"full\"` fits.")
  }
  if (length(SUBJECT_ID) != 1 || !SUBJECT_ID %in% FIT$data$labs$obs) {
    stop("`SUBJECT_ID` must match exactly one entry in `FIT$data$labs$obs`.")
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
  subject_index <- match(SUBJECT_ID, FIT$data$labs$obs)

  grades <- FIT$data$gradesMat
  times <- FIT$data$timeMat
  todo <- FIT$data$todoMat
  outcome_full <- FIT$data$outcome
  last_year_full <- FIT$data$last_year
  max_time_full <- FIT$data$max_time
  covariates <- as.matrix(FIT$data$X)

  subject_times <- times[subject_index, ]
  subject_todo <- todo[subject_index, ] == 1

  exam_days <- subject_times[subject_todo]
  exam_names <- FIT$data$labs$exams[subject_todo]
  order_idx <- order(exam_days, na.last = NA)
  event_days <- exam_days[order_idx]
  event_exams <- exam_names[order_idx]

  cutoff_days <- c(0, as.numeric(event_days))
  trigger_exam <- c(NA_character_, as.character(event_exams))

  theta <- FIT$fit$par
  dim_irt <- data_dims$n_exams * (data_dims$n_grades + 3)
  n_cov <- data_dims$n_cov

  if (is.null(horizons)) {
    horizons <- seq_len(yb)
  }

  compute_estep <- function(cutoff_day) {
    grades_trunc <- grades
    times_trunc <- times

    mask_future <- !is.na(times_trunc[subject_index, ]) &
      times_trunc[subject_index, ] > cutoff_day
    grades_trunc[subject_index, mask_future] <- NA_integer_
    times_trunc[subject_index, mask_future] <- NA

    max_time_trunc <- max_time_full
    exam_limit <- if (cutoff_day > 0) cutoff_day else 1L
    max_time_trunc[subject_index] <- as.integer(pmax(
      pmin(max_time_full[subject_index], exam_limit),
      1L
    ))

    cutoff_year_floor <- floor(cutoff_day / YEAR_LENGTH)
    last_year_trunc <- last_year_full
    outcome_trunc <- outcome_full
    last_year_trunc[subject_index] <- as.integer(pmin(
      last_year_full[subject_index],
      cutoff_year_floor
    ))
    outcome_trunc[subject_index] <- if (
      last_year_full[subject_index] > cutoff_year_floor
    ) {
      0L
    } else {
      outcome_full[subject_index]
    }

    yle_trunc <- FIT$data$yle
    subject_required <- todo[subject_index, ]
    obs_times <- times_trunc[subject_index, subject_required]
    if (length(obs_times) > 0 && any(!is.na(obs_times))) {
      yle_trunc[subject_index] <- as.integer(ceiling(
        max(obs_times, na.rm = TRUE) / YEAR_LENGTH
      ))
    } else {
      yle_trunc[subject_index] <- 100L
    }

    estep_out <- cpp_estep(
      THETA = theta,
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
    list(
      Ew = estep_out$Ew,
      nodes = estep_out$nodes,
      grades_trunc = grades_trunc,
      times_trunc = times_trunc,
      cutoff_year_floor = cutoff_year_floor
    )
  }

  summarize_event <- function(event_index, cutoff_day, trigger_exam) {
    estep <- compute_estep(cutoff_day)
    Ew <- estep$Ew
    nodes <- estep$nodes
    grades_trunc <- estep$grades_trunc
    times_trunc <- estep$times_trunc
    cutoff_year_floor <- estep$cutoff_year_floor

    Ew[!is.finite(Ew)] <- 0
    row_sums <- rowSums(Ew)
    zero_rows <- row_sums <= 0
    nq <- nrow(nodes)
    if (any(zero_rows)) {
      Ew[zero_rows, ] <- 1 / nq
      row_sums[zero_rows] <- 1
    }
    Ew <- Ew / row_sums

    mu_ability <- mu_speed <- rep(0, nrow(covariates))
    if (n_cov > 0) {
      gamma_ability <- theta[(dim_irt + 3):(dim_irt + 2 + n_cov)]
      gamma_speed <- theta[(dim_irt + 3 + n_cov):(dim_irt + 2 + 2 * n_cov)]
      mu_ability <- as.numeric(covariates %*% gamma_ability)
      mu_speed <- as.numeric(covariates %*% gamma_speed)
    }

    target_years <- horizons[horizons > cutoff_year_floor]
    if (length(target_years) == 0) {
      return(tibble())
    }
    target_days <- as.integer(round(
      (target_years - 1L) * YEAR_LENGTH + GRAD_EXTENSION * YEAR_LENGTH
    ))

    Ew_row <- Ew[subject_index, ]
    cov_i <- if (n_cov > 0) {
      as.numeric(covariates[subject_index, ])
    } else {
      numeric(0)
    }
    pending_idx <- which(
      todo[subject_index, ] & is.na(grades_trunc[subject_index, ])
    )
    finished <- length(pending_idx) == 0

    ability_eap <- sum(Ew_row * (nodes[, 1] + mu_ability[subject_index]))
    speed_eap <- sum(Ew_row * (nodes[, 2] + mu_speed[subject_index]))

    prob_grad_years <- numeric(length(target_years))
    prob_drop_years <- numeric(length(target_years))
    prob_trans_years <- numeric(length(target_years))

    for (q in seq_len(nq)) {
      w <- Ew_row[q]
      if (w <= 0 || !is.finite(w)) {
        next
      }
      ability <- nodes[q, 1] + mu_ability[subject_index]
      speed <- nodes[q, 2] + mu_speed[subject_index]

      hazard_grad_vec <- hazard_drop_vec <- hazard_trans_vec <- numeric(length(
        target_years
      ))
      for (t_idx in seq_along(target_years)) {
        ty <- target_years[t_idx]
        hazard_grad_vec[t_idx] <- cpp_hazard(
          1L,
          ty,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
        hazard_drop_vec[t_idx] <- cpp_hazard(
          2L,
          ty,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
        hazard_trans_vec[t_idx] <- cpp_hazard(
          3L,
          ty,
          theta,
          cov_i,
          ability,
          speed,
          data_dims$yb,
          FALSE
        )$prob
      }

      if (finished) {
        prob_grad_years <- prob_grad_years + w * hazard_grad_vec
        next
      }

      prob_finish_vec <- rep(1, length(target_years))
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
        for (t_idx in seq_along(target_years)) {
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
          if (S_target <= 1e-12) {
            S_target <- 1e-12
          }
          prob_exam <- 1 - (S_target / S_cut)
          prob_exam <- min(max(prob_exam, 0), 1)
          prob_finish_vec[t_idx] <- prob_finish_vec[t_idx] * prob_exam
        }
      }

      prob_finish_vec <- pmin(pmax(prob_finish_vec, 0), 1)
      prob_grad_years <- prob_grad_years +
        w * (prob_finish_vec * hazard_grad_vec)
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
    prob_grad_cum <- numeric(length(target_years))
    prob_drop_cum <- numeric(length(target_years))
    prob_trans_cum <- numeric(length(target_years))
    prob_still <- numeric(length(target_years))

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

    occurred_before_cutoff <- outcome_full %in%
      c(1L, 2L, 3L) &
      last_year_full <= cutoff_year_floor
    if (occurred_before_cutoff[subject_index]) {
      if (outcome_full[subject_index] == 1L) {
        prob_grad_cum <- rep(1, length(target_years))
        prob_drop_cum <- prob_trans_cum <- prob_still <- rep(
          0,
          length(target_years)
        )
      } else if (outcome_full[subject_index] == 2L) {
        prob_drop_cum <- rep(1, length(target_years))
        prob_grad_cum <- prob_trans_cum <- prob_still <- rep(
          0,
          length(target_years)
        )
      } else if (outcome_full[subject_index] == 3L) {
        prob_trans_cum <- rep(1, length(target_years))
        prob_grad_cum <- prob_drop_cum <- prob_still <- rep(
          0,
          length(target_years)
        )
      }
    }

    tibble(
      event_index = event_index,
      trigger_exam = trigger_exam,
      trigger_day = cutoff_day,
      cutoff_day = cutoff_day,
      cutoff_year = cutoff_year_floor,
      target_year = target_years,
      prob_graduation = prob_grad_cum,
      prob_dropout = prob_drop_cum,
      prob_transfer = prob_trans_cum,
      prob_still_enrolled = prob_still,
      ability_eap = ability_eap,
      speed_eap = speed_eap,
      n_exams_done = sum(!is.na(times_trunc[subject_index, ])),
      n_exams_todo = sum(todo[subject_index, ])
    )
  }

  purrr::map_dfr(
    seq_along(cutoff_days),
    ~ summarize_event(.x - 1L, cutoff_days[.x], trigger_exam[.x])
  ) |>
    arrange(.data$event_index, .data$target_year)
}
