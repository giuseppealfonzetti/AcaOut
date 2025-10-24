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
#'
#' @return A list structured as the output of [check_data()], enriched with the
#'   true parameter list (element `params`) and the generated latent scores
#'   (element `latent`). The object can be passed directly to [fit_EM()] or
#'   [fit_BFGS()].
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats plogis qlogis runif rnorm
simulate_crirt_data <- function(
  N_STUDENTS = 100,
  N_EXAMS = 10,
  N_GRADES = 4,
  MAX_YEAR = 5,
  N_COV = 2,
  SEED = 123,
  PARAMS = NULL,
  LATMAT = NULL
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
      irt[j, ] <- c(thresholds, slope, zeta, omega)
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

  stopifnot(nrow(params$irt) == N_EXAMS)
  stopifnot(ncol(params$irt) == N_GRADES + 3)
  stopifnot(length(params$cr$grad) == 1)
  stopifnot(nrow(params$cr$beta) == MAX_YEAR + 2)
  stopifnot(ncol(params$cr$beta) == 2)

  if (is.null(params$lat_reg)) {
    params$lat_reg <- matrix(0, nrow = 2, ncol = N_COV)
  } else {
    stopifnot(ncol(params$lat_reg) == N_COV)
  }

  covariates <- matrix(rnorm(N_STUDENTS * N_COV), ncol = N_COV)
  colnames(covariates) <- paste0("cov_", seq_len(N_COV))
  if (N_COV >= 1) {
    covariates[, 1] <- rbinom(N_STUDENTS, 1, 0.5)
  }

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
  todo_mat <- matrix(TRUE, nrow = N_STUDENTS, ncol = N_EXAMS)
  colnames(grades_mat) <- colnames(time_mat) <- colnames(todo_mat) <- paste0(
    "exam_",
    seq_len(N_EXAMS)
  )
  rownames(grades_mat) <- rownames(time_mat) <- rownames(todo_mat) <- paste0(
    "student_",
    seq_len(N_STUDENTS)
  )

  max_day <- MAX_YEAR * 365L
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

    # Potential exam outcomes (before censoring)
    potential_grades <- integer(N_EXAMS)
    potential_times <- numeric(N_EXAMS)

    for (j in seq_len(N_EXAMS)) {
      irt_row <- params$irt[j, ]
      thresholds <- irt_row[seq_len(N_GRADES)]
      slope <- irt_row[N_GRADES + 1]
      zeta <- irt_row[N_GRADES + 2]
      omega <- irt_row[N_GRADES + 3]

      logits_ge <- plogis(slope * (ability - thresholds))
      probs <- numeric(N_GRADES)
      for (k in seq_len(N_GRADES)) {
        upper <- if (k == N_GRADES) {
          0
        } else {
          plogis(slope * (ability - thresholds[k + 1]))
        }
        probs[k] <- logits_ge[k] - upper
      }
      probs <- pmax(probs, 1e-8)
      probs <- probs / sum(probs)
      potential_grades[j] <- sample.int(N_GRADES, size = 1, prob = probs)

      mu_time <- zeta - speed
      sd_time <- 1 / max(omega, 1e-3)
      potential_times[j] <- round(exp(rnorm(1, mean = mu_time, sd = sd_time)))
      potential_times[j] <- max(1, potential_times[j])
    }

    complete_day <- max(potential_times)
    year_complete <- if (complete_day <= max_day) {
      ceiling(complete_day / 365)
    } else {
      Inf
    }

    draw_event <- function() {
      latent_vec <- c(ability, speed)
      for (yr in seq_len(MAX_YEAR)) {
        if (yr < year_complete || is.infinite(year_complete)) {
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
          p_g <- plogis(params$cr$grad)
          if (runif(1) < p_g) {
            return(list(type = 1L, year = yr))
          }
        }
      }
      list(type = 0L, year = MAX_YEAR)
    }

    event <- draw_event()
    outcomes[i] <- event$type
    last_year[i] <- if (event$type == 0L) MAX_YEAR else event$year

    if (event$type == 1L) {
      event_day <- complete_day
      observed_times <- potential_times
      last_exam_year[i] <- as.integer(ceiling(complete_day / 365))
    } else if (event$type == 0L) {
      event_day <- max_day
      observed_times <- potential_times
      observed_times[potential_times > event_day] <- NA
    } else {
      event_day <- event$year * 365
      observed_times <- potential_times
      observed_times[potential_times > event_day] <- NA
    }

    observed_grades <- potential_grades
    observed_grades[is.na(observed_times)] <- NA

    grades_mat[i, ] <- observed_grades
    time_mat[i, ] <- observed_times
    max_time[i] <- as.integer(min(event_day, max_day))

    if (all(is.na(observed_times))) {
      last_exam_year[i] <- NA_integer_
    } else {
      last_exam_year[i] <- as.integer(ceiling(
        max(observed_times, na.rm = TRUE) / 365
      ))
    }
  }

  first_year <- rep(1L, N_STUDENTS)
  student_ids <- rownames(grades_mat)
  exam_ids <- colnames(grades_mat)
  grade_ids <- paste0("grade_", seq_len(N_GRADES))

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
