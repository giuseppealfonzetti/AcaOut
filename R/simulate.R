#' Simulate student-level data for CR-IRT models
#'
#' @param N_STUDENTS Number of students to generate.
#' @param N_EXAMS Number of exams in the study plan.
#' @param N_GRADES Number of ordered grade categories (minimum 2).
#' @param MAX_YEAR Maximum observed academic year (competing risks apply up to this year).
#' @param N_COV Number of student-level covariates.
#' @param SEED Random seed used for reproducibility.
#' @param PARAMS Optional list of model parameters as returned by [parVec2List()].
#' @param LATMAT Optional `N_STUDENTS x 2` matrix of latent ability and speed scores.
#' @param TODO Optional `N_STUDENTS x N_EXAMS` study plan; exams left out stay missing.
#' @param X Optional `N_STUDENTS x N_COV` covariate matrix. Generated when missing.
#' @param FIRST_YEAR Optional enrolment years. Dropout and transfer start from that year.
#' @param ADMIN_YEAR Optional last year of follow up per student. Defaults to `MAX_YEAR`.
#' @param YEAR_COMPLETE Optional year each student completes the plan, `100` when never.
#'   Derived from the simulated exam times when missing.
#' @param MAX_TIME_OFFSET Optional days added to `last_year * 365` for the last observable
#'   day. Defaults to zero.
#' @param YEAR_OFFSET Day of the year separating two academic years. Defaults to zero.
#' @param CHECK Set to FALSE to skip [check_data()] and return the raw simulated bundle.
#'
#' @return A list structured as the output of [check_data()], with the
#'   true parameter list (element `params`) and the generated latent scores
#'   (element `latent`). The object can be passed directly to [fit_EM()] or
#'   [fit_BFGS()].
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats plogis qlogis runif rnorm rbinom
simulate_crirt_data <- function(
  N_STUDENTS = 100,
  N_EXAMS = 10,
  N_GRADES = 4,
  MAX_YEAR = 5,
  N_COV = 2,
  SEED = 123,
  PARAMS = NULL,
  LATMAT = NULL,
  TODO = NULL,
  X = NULL,
  FIRST_YEAR = NULL,
  ADMIN_YEAR = NULL,
  YEAR_COMPLETE = NULL,
  MAX_TIME_OFFSET = NULL,
  YEAR_OFFSET = 0L,
  CHECK = TRUE
) {
  stopifnot(N_STUDENTS > 0)
  stopifnot(N_EXAMS > 0)
  stopifnot(N_GRADES > 1)
  stopifnot(MAX_YEAR > 0)
  stopifnot(N_COV >= 1)

  set.seed(SEED)

  log_competing <- function(p, other) {
    rest <- 1 - p - other
    if (rest <= 0) {
      rest <- 1e-3
    }
    log(p / rest)
  }

  generate_default_params <- function() {
    irt <- matrix(0, nrow = N_EXAMS, ncol = N_GRADES + 3)
    for (j in seq_len(N_EXAMS)) {
      base_threshold <- seq(-1.8, 1.5, length.out = N_GRADES)
      thresholds <- sort(base_threshold + rnorm(N_GRADES, 0, 0.25))
      slope <- runif(1, 0.8, 5)
      base_time <- log(130 + 25 * j)
      zeta <- rnorm(1, base_time, 0.12)
      omega <- runif(1, 0.9, 1.4)
      irt[j, ] <- c(slope * thresholds, slope, zeta, omega)
    }

    speed_sd <- runif(1, 0.55, 0.75)
    rho <- runif(1, 0.15, 0.35)
    lat_var <- matrix(
      c(1, rho * speed_sd, rho * speed_sd, speed_sd^2),
      nrow = 2
    )

    lat_reg <- matrix(0, nrow = 2, ncol = N_COV)
    lat_reg[1, ] <- rnorm(N_COV, 0, 0.25)
    lat_reg[2, ] <- rnorm(N_COV, -0.1, 0.2)

    dropout_target <- seq(0.18, 0.05, length.out = MAX_YEAR)
    transfer_target <- seq(0.08, 0.02, length.out = MAX_YEAR)
    beta_mat <- matrix(0, nrow = MAX_YEAR + 2, ncol = 2)
    for (yr in seq_len(MAX_YEAR)) {
      beta_mat[yr, 1] <- log_competing(dropout_target[yr], transfer_target[yr])
      beta_mat[yr, 2] <- log_competing(transfer_target[yr], dropout_target[yr])
    }
    beta_mat[MAX_YEAR + 1, ] <- c(-0.8, -0.5) # ability effects
    beta_mat[MAX_YEAR + 2, ] <- c(-0.7, -0.6) # speed effects

    grad_intercept <- qlogis(0.62)

    list(
      irt = irt,
      lat_var = lat_var,
      lat_reg = if (N_COV > 0) lat_reg else NULL,
      cr = list(
        grad = grad_intercept,
        beta = beta_mat
      )
    )
  }

  params <- PARAMS
  if (is.null(params)) {
    params <- generate_default_params()
  }

  grad_intercept <- if (!is.null(params$cr$grad)) {
    params$cr$grad
  } else {
    params$cr$graduation
  }

  stopifnot(nrow(params$irt) == N_EXAMS)
  stopifnot(ncol(params$irt) == N_GRADES + 3)
  stopifnot(length(grad_intercept) == 1)
  stopifnot(nrow(params$cr$beta) == MAX_YEAR + 2)
  stopifnot(ncol(params$cr$beta) == 2)

  if (is.null(params$lat_reg)) {
    params$lat_reg <- matrix(0, nrow = 2, ncol = N_COV)
  } else {
    stopifnot(ncol(params$lat_reg) == N_COV)
  }

  covariates <- if (is.null(X)) {
    cov_mat <- matrix(rnorm(N_STUDENTS * N_COV), ncol = N_COV)
    if (N_COV >= 1) {
      cov_mat[, 1] <- rbinom(N_STUDENTS, 1, 0.5)
    }
    colnames(cov_mat) <- paste0("cov_", seq_len(N_COV))
    cov_mat
  } else {
    cov_mat <- as.matrix(X)
    stopifnot(nrow(cov_mat) == N_STUDENTS, ncol(cov_mat) == N_COV)
    cov_mat
  }

  todo_mat <- if (is.null(TODO)) {
    matrix(TRUE, nrow = N_STUDENTS, ncol = N_EXAMS)
  } else {
    matrix(as.logical(TODO), nrow = N_STUDENTS, ncol = N_EXAMS)
  }

  first_year <- if (is.null(FIRST_YEAR)) {
    rep(1L, N_STUDENTS)
  } else {
    as.integer(FIRST_YEAR)
  }

  admin_year <- if (is.null(ADMIN_YEAR)) {
    rep(as.integer(MAX_YEAR), N_STUDENTS)
  } else {
    as.integer(ADMIN_YEAR)
  }

  time_offset <- if (is.null(MAX_TIME_OFFSET)) {
    rep(0L, N_STUDENTS)
  } else {
    as.integer(MAX_TIME_OFFSET)
  }

  stopifnot(length(first_year) == N_STUDENTS)
  stopifnot(length(admin_year) == N_STUDENTS)
  stopifnot(length(time_offset) == N_STUDENTS)
  stopifnot(all(admin_year <= MAX_YEAR))

  year_of <- function(t) pmax(1L, as.integer(ceiling((t - YEAR_OFFSET) / 365)))

  latent_mat <- if (is.null(LATMAT)) {
    mu_mat <- covariates %*% t(params$lat_reg)
    latent_draws <- matrix(NA_real_, nrow = N_STUDENTS, ncol = 2)
    for (i in seq_len(N_STUDENTS)) {
      latent_draws[i, ] <- MASS::mvrnorm(
        1,
        mu = c(mu_mat[i, ]),
        Sigma = params$lat_var
      )
    }
    latent_draws
  } else {
    stopifnot(nrow(LATMAT) == N_STUDENTS, ncol(LATMAT) == 2)
    LATMAT
  }
  colnames(latent_mat) <- c("ability", "speed")

  grades_mat <- matrix(NA_integer_, nrow = N_STUDENTS, ncol = N_EXAMS)
  time_mat <- matrix(NA_integer_, nrow = N_STUDENTS, ncol = N_EXAMS)
  colnames(grades_mat) <- colnames(time_mat) <- colnames(todo_mat) <- paste0(
    "exam_",
    seq_len(N_EXAMS)
  )
  rownames(grades_mat) <- rownames(time_mat) <- rownames(todo_mat) <- paste0(
    "student_",
    seq_len(N_STUDENTS)
  )

  outcomes <- integer(N_STUDENTS)
  last_year <- rep(1L, N_STUDENTS)
  last_exam_year <- rep(100L, N_STUDENTS)
  max_time <- integer(N_STUDENTS)

  ability_effects <- params$cr$beta[MAX_YEAR + 1, ]
  speed_effects <- params$cr$beta[MAX_YEAR + 2, ]
  year_intercepts <- params$cr$beta[seq_len(MAX_YEAR), , drop = FALSE]

  for (i in seq_len(N_STUDENTS)) {
    ability <- latent_mat[i, 1]
    speed <- latent_mat[i, 2]
    plan <- which(todo_mat[i, ])

    potential_grades <- rep(NA_integer_, N_EXAMS)
    potential_times <- rep(NA_real_, N_EXAMS)

    for (j in plan) {
      irt_row <- params$irt[j, ]
      intercepts <- irt_row[seq_len(N_GRADES)]
      slope <- irt_row[N_GRADES + 1]
      zeta <- irt_row[N_GRADES + 2]
      omega <- irt_row[N_GRADES + 3]

      cum <- plogis(slope * ability - intercepts)
      if (runif(1) >= cum[1]) {
        potential_times[j] <- Inf
        next
      }

      probs <- cum - c(cum[-1], 0)
      probs <- pmax(probs, 1e-8)
      probs <- probs / sum(probs)
      potential_grades[j] <- sample.int(N_GRADES, size = 1, prob = probs)

      mu_time <- zeta - speed
      sd_time <- 1 / max(omega, 1e-3)
      potential_times[j] <- max(1, round(exp(rnorm(1, mean = mu_time, sd = sd_time))))
    }

    year_complete <- if (!is.null(YEAR_COMPLETE)) {
      if (YEAR_COMPLETE[i] >= 100) Inf else as.integer(YEAR_COMPLETE[i])
    } else {
      complete_day <- if (length(plan) > 0) max(potential_times[plan]) else NA_real_
      if (isTRUE(is.finite(complete_day))) year_of(complete_day) else Inf
    }
    if (is.finite(year_complete) && year_complete > admin_year[i]) {
      year_complete <- Inf
    }

    draw_event <- function() {
      if (first_year[i] > admin_year[i]) {
        return(list(type = 0L, year = admin_year[i]))
      }
      for (yr in seq.int(first_year[i], admin_year[i])) {
        if (yr < year_complete) {
          eta_d <- year_intercepts[yr, 1] +
            ability_effects[1] * ability +
            speed_effects[1] * speed
          eta_t <- year_intercepts[yr, 2] +
            ability_effects[2] * ability +
            speed_effects[2] * speed
          den <- 1 + exp(eta_d) + exp(eta_t)
          p_d <- exp(eta_d) / den
          p_t <- exp(eta_t) / den
          draw <- runif(1)
          if (draw < p_d) {
            return(list(type = 2L, year = yr))
          }
          if (draw < p_d + p_t) {
            return(list(type = 3L, year = yr))
          }
        } else {
          if (runif(1) < plogis(grad_intercept)) {
            return(list(type = 1L, year = yr))
          }
        }
      }
      list(type = 0L, year = admin_year[i])
    }

    event <- draw_event()
    outcomes[i] <- event$type
    last_year[i] <- if (event$type == 0L) admin_year[i] else event$year

    censor_day <- last_year[i] * 365 + time_offset[i]
    observed_times <- potential_times
    observed_times[!is.na(potential_times) & potential_times > censor_day] <- NA

    observed_grades <- potential_grades
    observed_grades[is.na(observed_times)] <- NA

    grades_mat[i, ] <- observed_grades
    time_mat[i, ] <- as.integer(observed_times)
    max_time[i] <- as.integer(censor_day)

    leaver <- event$type %in% c(2L, 3L)
    last_exam_year[i] <- if (
      !is.finite(year_complete) || (is.null(YEAR_COMPLETE) && leaver)
    ) {
      100L
    } else {
      as.integer(year_complete)
    }
  }

  student_ids <- rownames(grades_mat)
  exam_ids <- colnames(grades_mat)
  grade_ids <- paste0("grade_", seq_len(N_GRADES))

  if (!CHECK) {
    return(list(
      gradesMat = grades_mat,
      timeMat = time_mat,
      todoMat = todo_mat,
      outcome = outcomes,
      first_year = first_year,
      last_year = last_year,
      yle = last_exam_year,
      max_time = max_time,
      X = covariates,
      params = params,
      latent = latent_mat
    ))
  }

  data_obj <- check_data(
    GRADES = grades_mat,
    TIMES = time_mat,
    TODO = todo_mat,
    OUTCOME = outcomes,
    X = covariates,
    FIRST_YEAR = first_year,
    LAST_YEAR = last_year,
    LAST_EXAM_YEAR = last_exam_year,
    MAX_TIME = max_time,
    LABS_EXAMS = exam_ids,
    LABS_OBS = student_ids,
    LABS_GRADES = grade_ids,
    LABS_COV = colnames(covariates),
    VERBOSE = FALSE
  )

  data_obj$n_exams <- N_EXAMS
  data_obj$n_grades <- N_GRADES
  data_obj$n_cov <- N_COV
  data_obj$yb <- MAX_YEAR

  data_obj$params <- params
  data_obj$latent <- latent_mat

  data_obj
}

#' Simulate a dataset from a fitted model, reusing its design
#'
#' Generates grades, times, outcomes and censoring from the parameters of `FIT`,
#' keeping the study plan, covariates, enrolment years and observation windows of
#' its students. Parameters are held at their estimates; nothing is re-estimated.
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS].
#' @param SEED Random seed used for reproducibility.
#' @param ADMIN_YEAR Last year of follow up per student. Defaults to the model's `yb`.
#' @param YEAR_OFFSET Day of the year separating two academic years.
#'
#' @return A list shaped like `FIT`, whose `data` element carries the simulated
#'   grades, times, outcomes, last years and completion years, ready for [compute_map()].
#'
#' @export
simulate_from_fit <- function(FIT, SEED, ADMIN_YEAR = NULL, YEAR_OFFSET = 0L) {
  np <- FIT$data$data_dims

  sim <- simulate_crirt_data(
    N_STUDENTS = np$n_obs,
    N_EXAMS = np$n_exams,
    N_GRADES = np$n_grades,
    MAX_YEAR = np$yb,
    N_COV = np$n_cov,
    SEED = SEED,
    PARAMS = parVec2List(
      FIT$fit$par,
      N_GRADES = np$n_grades,
      N_EXAMS = np$n_exams,
      N_COV = np$n_cov,
      YB = np$yb
    ),
    TODO = FIT$data$todoMat,
    X = FIT$data$X,
    FIRST_YEAR = FIT$data$first_year,
    ADMIN_YEAR = ADMIN_YEAR,
    MAX_TIME_OFFSET = FIT$data$max_time - FIT$data$last_year * 365,
    YEAR_OFFSET = YEAR_OFFSET,
    CHECK = FALSE
  )

  dat <- FIT$data
  dimnames(sim$gradesMat) <- dimnames(dat$gradesMat)
  dimnames(sim$timeMat) <- dimnames(dat$timeMat)
  dat$gradesMat <- sim$gradesMat
  dat$timeMat <- sim$timeMat
  dat$outcome <- stats::setNames(sim$outcome, names(dat$outcome))
  dat$last_year <- stats::setNames(sim$last_year, names(dat$last_year))
  dat$yle <- stats::setNames(sim$yle, names(dat$yle))
  dat$max_time <- stats::setNames(sim$max_time, names(dat$max_time))

  out <- FIT
  out$data <- dat
  out$latent <- sim$latent
  out
}
