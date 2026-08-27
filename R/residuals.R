#' Grade and time residuals conditioned on plug-in latent traits
#'
#' Residuals of the measurement model at a single value of the latent traits
#' rather than integrating them out. Each residual is the observation minus its
#' conditional expectation over the set where the exam was observed; numerators
#' and standard deviations are returned separately, so the standardised residual
#' is `grade / grade_sd` and `time / time_sd`. With the MAP scores as plug-in,
#' the cross-exam correlations of the grade residuals are Yen's Q3 statistics.
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS].
#' @param LATMAT Optional `n_obs x 2` matrix of ability and speed to plug in.
#'   Defaults to the MAP scores from [compute_map()].
#' @param LAST_DAY Optional last observed day per student, at which the times are
#'   truncated. Defaults to `FIT$data$max_time`.
#'
#' @return A list of four `n_obs x n_exams` matrices: `grade` and `time` hold the
#'   observation minus its conditional expectation, `grade_sd` and `time_sd` the
#'   conditional standard deviations. `grade` and `time` follow the observation
#'   pattern of the data; the standard deviations are given wherever they exist.
#'
#' @export
conditioned_residuals <- function(FIT, LATMAT = NULL, LAST_DAY = NULL) {
  np <- FIT$data$data_dims
  pars <- parVec2List(
    FIT$fit$par,
    N_GRADES = np$n_grades,
    N_EXAMS = np$n_exams,
    N_COV = np$n_cov,
    YB = np$yb
  )

  lat <- if (is.null(LATMAT)) {
    as.matrix(suppressMessages(compute_map(FIT, TIDY = FALSE)))
  } else {
    LATMAT
  }

  last_day <- if (is.null(LAST_DAY)) FIT$data$max_time else LAST_DAY
  stopifnot(length(last_day) == np$n_obs)
  log_c <- log(last_day)

  grade_res <- time_res <- grade_sd <- time_sd <- matrix(
    NA_real_,
    np$n_obs,
    np$n_exams,
    dimnames = list(NULL, FIT$data$labs$exams)
  )

  for (j in seq_len(np$n_exams)) {
    cum <- vapply(
      seq_len(np$n_grades),
      function(g) stats::plogis(pars$irt[j, "slope"] * lat[, 1] - pars$irt[j, g]),
      numeric(np$n_obs)
    )
    m1 <- rowSums(cum) / cum[, 1]
    m2 <- drop(cum %*% (2 * seq_len(np$n_grades) - 1)) / cum[, 1]
    v_grade <- m2 - m1^2

    mu <- pars$irt[j, "time_loc"] - lat[, 2]
    om <- pars$irt[j, "time_invsd"]
    kap <- om * (log_c - mu)
    lam <- exp(stats::dnorm(kap, log = TRUE) - stats::pnorm(kap, log.p = TRUE))
    v_time <- (1 - kap * lam - lam^2) / om^2

    v_grade[!(v_grade > 0)] <- NA_real_
    v_time[!(v_time > 0)] <- NA_real_

    grade_res[, j] <- FIT$data$gradesMat[, j] - m1
    time_res[, j] <- log(FIT$data$timeMat[, j]) - (mu - lam / om)
    grade_sd[, j] <- sqrt(v_grade)
    time_sd[, j] <- sqrt(v_time)
  }

  out <- list(
    grade = grade_res,
    time = time_res,
    grade_sd = grade_sd,
    time_sd = time_sd
  )
  lapply(out, function(m) {
    m[!is.finite(m)] <- NA
    m
  })
}

#' Residuals of datasets simulated under local independence
#'
#' Simulates `B` datasets from `FIT` with [simulate_from_fit()], which satisfy
#' local independence by construction, re-scores the latent traits on each and
#' returns the residuals. Nothing is summarised: correlation panels and tests are
#' left to the caller. By default parameters are held at `FIT$fit$par` and only
#' the latent scores are re-estimated; supplying `REFIT` re-estimates the
#' parameters too, warm starting [fit_BFGS()] at `FIT$fit$par`.
#'
#' The returned object is large: with `n_obs = 786`, `n_exams = 22` and
#' `B = 1000` the four arrays occupy about 553 MB in memory.
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS].
#' @param B Number of replications.
#' @param NCORES Number of cores. Defaults to `detectCores() - 1`, or `1` on Windows.
#' @param SEED Replication `b` is simulated with seed `SEED + b`.
#' @param ADMIN_YEAR Last year of follow up per student, passed to [simulate_from_fit()].
#' @param YEAR_OFFSET Day of the year separating two academic years.
#' @param REFIT Optional control list for [ucminf::ucminf()], passed to [fit_BFGS()].
#'   When given, each replication is refitted from `FIT$fit$par`; the convergence
#'   codes are returned in `meta`.
#' @param VERBOSE TRUE to report progress.
#'
#' @return A list with `observed` (residuals of the data `FIT` was estimated on),
#'   `simulated` (the `n_obs x n_exams x B` arrays `grade`, `time`, `grade_sd`,
#'   `time_sd`) and `meta` (the generating arguments and, with `REFIT`, the
#'   convergence code of each replication).
#'
#' @export
simulate_residuals <- function(
  FIT,
  B = 1000,
  NCORES = NULL,
  SEED = 1,
  ADMIN_YEAR = NULL,
  YEAR_OFFSET = 0L,
  REFIT = NULL,
  VERBOSE = FALSE
) {
  ncores <- if (is.null(NCORES)) {
    if (.Platform$OS.type == "windows") {
      1L
    } else {
      max(1L, parallel::detectCores() - 1L)
    }
  } else {
    as.integer(NCORES)
  }

  np <- FIT$data$data_dims
  exams <- FIT$data$labs$exams

  one_rep <- function(b) {
    sim <- simulate_from_fit(
      FIT,
      SEED = SEED + b,
      ADMIN_YEAR = ADMIN_YEAR,
      YEAR_OFFSET = YEAR_OFFSET
    )

    conv <- NA_integer_
    if (!is.null(REFIT)) {
      utils::capture.output(
        refit <- suppressMessages(fit_BFGS(
          DATA = sim$data,
          GRID = FIT$grid,
          WEIGHTS = FIT$weights,
          THETA_START = FIT$fit$par,
          MOD = FIT$mod,
          control = REFIT
        )),
        file = nullfile()
      )
      conv <- refit$fit$convergence
      if (conv < 0) {
        stop(
          "replication ", b, ": ucminf did not start (convergence ", conv, "): ",
          refit$fit$message
        )
      }
      sim$fit$par <- refit$fit$par
    }

    res <- conditioned_residuals(sim)
    attr(res, "convergence") <- conv
    res
  }

  if (ncores > 1L && interactive()) {
    warning(
      "Forking is unreliable inside a graphical R session. Run this from Rscript, ",
      "or set NCORES = 1.",
      call. = FALSE
    )
  }

  if (VERBOSE) {
    message("Simulating ", B, " datasets on ", ncores, " cores.")
  }

  reps <- if (ncores > 1L) {
    parallel::mclapply(seq_len(B), one_rep, mc.cores = ncores)
  } else {
    lapply(seq_len(B), one_rep)
  }

  died <- vapply(reps, is.null, logical(1))
  errored <- vapply(reps, function(x) inherits(x, "try-error"), logical(1))
  if (any(died) || any(errored)) {
    reason <- if (any(errored)) {
      conditionMessage(attr(reps[[which(errored)[1]]], "condition"))
    } else {
      paste(
        "the workers died without returning a result, which usually means the",
        "replications were forked from a graphical R session. Run this from Rscript,",
        "or set NCORES = 1."
      )
    }
    stop(
      sum(died | errored), " of ", B, " replications failed. ", reason,
      call. = FALSE
    )
  }

  template <- matrix(0, np$n_obs, np$n_exams)
  families <- c("grade", "time", "grade_sd", "time_sd")
  simulated <- lapply(stats::setNames(families, families), function(what) {
    arr <- vapply(reps, function(x) x[[what]], template)
    dimnames(arr) <- list(NULL, exams, NULL)
    arr
  })

  list(
    observed = conditioned_residuals(FIT),
    simulated = simulated,
    meta = list(
      B = B,
      seed = SEED,
      admin_year = ADMIN_YEAR,
      year_offset = YEAR_OFFSET,
      exams = exams,
      refit = REFIT,
      convergence = vapply(reps, function(x) attr(x, "convergence"), integer(1))
    )
  )
}
