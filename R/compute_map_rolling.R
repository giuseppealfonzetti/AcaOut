#' Rolling MAP latent scores
#'
#' Compute MAP estimates of ability and speed after truncating the history at a
#' given academic year. The truncation matches the logic used by
#' [predict_one_year()], i.e. grades and times after
#' `CUTOFF_YEAR * YEAR_LENGTH + (GRAD_EXTENSION - 1) * YEAR_LENGTH` are discarded and outcomes
#' occurring after `CUTOFF_YEAR` are recoded as "still enrolled".
#'
#' @inheritParams compute_map
#' @param CUTOFF_YEAR Academic year `y` whose information is used for the MAP
#'   calculation. Must satisfy `0 <= CUTOFF_YEAR < FIT$data$data_dims$yb`.
#' @param YEAR_LENGTH Length (in days) of an academic year. Defaults to 365.
#' @param GRAD_EXTENSION One plus the extra number of academic years retained
#'   beyond `CUTOFF_YEAR` when censoring exam information. Defaults to 1.5,
#'   meaning half an extra year is kept.
#' @param TIDY Logical; return a tibble instead of a matrix.
#'
#' @return A matrix (or tibble) with MAP estimates of ability and speed based on
#'   the truncated data.
#' @export
#'
#' @importFrom dplyr as_tibble
compute_map_rolling <- function(
  FIT,
  CUTOFF_YEAR,
  MATSTART = NULL,
  YEAR_LENGTH = 365,
  GRAD_EXTENSION = 1.5,
  TIDY = TRUE
) {
  if (!(FIT$mod %in% c("full", "grtc"))) {
    stop("Model not available. Provide fit object for `full` or `grtc` models.")
  }
  if (length(CUTOFF_YEAR) != 1 || !is.finite(CUTOFF_YEAR)) {
    stop("`CUTOFF_YEAR` must be a single finite value.")
  }
  CUTOFF_YEAR <- as.integer(CUTOFF_YEAR)
  if (CUTOFF_YEAR < 0) {
    stop("`CUTOFF_YEAR` must be non-negative.")
  }
  yb <- FIT$data$data_dims$yb
  if (CUTOFF_YEAR > yb) {
    stop(
      "`CUTOFF_YEAR` must be less than or equal to the maximum observed year."
    )
  }
  # if (CUTOFF_YEAR == yb) {
  #   return(compute_map(FIT, MATSTART = MATSTART, TIDY = TIDY))
  # }
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

  grades <- FIT$data$gradesMat
  times <- FIT$data$timeMat
  todo <- FIT$data$todoMat
  outcome_full <- FIT$data$outcome
  last_year_full <- FIT$data$last_year
  max_time_full <- FIT$data$max_time
  covariates <- as.matrix(FIT$data$X)
  first_year <- FIT$data$first_year
  labs <- FIT$data$labs
  n_students <- nrow(grades)

  cutoff_exam_day <- if (CUTOFF_YEAR == 0L) {
    0L
  } else {
    extra <- (GRAD_EXTENSION - 1) * YEAR_LENGTH
    as.integer(round(CUTOFF_YEAR * YEAR_LENGTH + extra))
  }

  grades_trunc <- grades
  times_trunc <- times
  mask_future <- !is.na(times_trunc) & times_trunc > cutoff_exam_day
  grades_trunc[mask_future] <- NA_integer_
  times_trunc[mask_future] <- NA_integer_

  exam_limit <- if (cutoff_exam_day > 0L) cutoff_exam_day else 1L
  max_time_trunc <- as.integer(pmax(pmin(max_time_full, exam_limit), 1L))
  last_year_trunc <- as.integer(pmin(last_year_full, CUTOFF_YEAR))

  outcome_trunc <- outcome_full
  outcome_trunc[last_year_full > CUTOFF_YEAR] <- 0L
  outcome_trunc <- as.integer(outcome_trunc)

  yle_trunc <- rep.int(100L, n_students)
  for (i in seq_len(n_students)) {
    required <- todo[i, ]
    obs_times <- times_trunc[i, required]
    if (length(obs_times) > 0 && all(!is.na(obs_times))) {
      yle_trunc[i] <- as.integer(ceiling(max(obs_times) / YEAR_LENGTH))
    }
  }

  if (is.null(MATSTART)) {
    MATSTART <- matrix(0, n_students, 2)
  } else {
    MATSTART <- as.matrix(MATSTART)
    if (!all(dim(MATSTART) == c(n_students, 2))) {
      stop(
        "`MATSTART` must be a matrix with `nrow(FIT$data$gradesMat)` rows and 2 columns."
      )
    }
  }

  if (FIT$mod == "full") {
    ID_COMPLETENLL <- function(LAT, ID) {
      -CRGRTCM_complete_obs(
        THETA = FIT$fit$par,
        EXAMS_GRADES = grades_trunc[ID, ],
        EXAMS_DAYS = times_trunc[ID, ],
        EXAMS_SET = todo[ID, ],
        EXAMS_OBSFLAG = !is.na(times_trunc[ID, ]),
        OUTCOME = outcome_trunc[ID],
        YEAR_FIRST = first_year[ID],
        YEAR_LAST = last_year_trunc[ID],
        YEAR_LAST_EXAM = yle_trunc[ID],
        COVARIATES = covariates[ID, ],
        MAX_DAY = max_time_trunc[ID],
        YB = yb,
        N_GRADES = FIT$data$data_dims$n_grades,
        N_EXAMS = FIT$data$data_dims$n_exams,
        ABILITY = LAT[1],
        SPEED = LAT[2]
      )
    }
  } else {
    ID_COMPLETENLL <- function(LAT, ID) {
      -GRTCM_complete_obs(
        THETA = FIT$fit$par,
        EXAMS_GRADES = grades_trunc[ID, ],
        EXAMS_DAYS = times_trunc[ID, ],
        EXAMS_SET = todo[ID, ],
        EXAMS_OBSFLAG = !is.na(times_trunc[ID, ]),
        COVARIATES = covariates[ID, ],
        MAX_DAY = max_time_trunc[ID],
        N_GRADES = FIT$data$data_dims$n_grades,
        N_EXAMS = FIT$data$data_dims$n_exams,
        ABILITY = LAT[1],
        SPEED = LAT[2]
      )
    }
  }

  map_mat <- Reduce(
    rbind,
    lapply(
      seq_len(n_students),
      function(ID) {
        fit_res <- ucminf::ucminf(
          par = as.numeric(MATSTART[ID, ]),
          fn = ID_COMPLETENLL,
          ID = ID
        )
        as.numeric(fit_res$par)
      }
    )
  )

  rownames(map_mat) <- labs$obs
  colnames(map_mat) <- c("ability", "speed")

  if (TIDY) {
    map_mat <- as_tibble(map_mat, rownames = "subject_id") |>
      mutate(year = CUTOFF_YEAR)
  }

  map_mat
}
