# synthetic fit helper
make_fit <- function(seed = 11, n_students = 60, params = NULL, check = TRUE) {
  n_exams <- 5L
  n_grades <- 4L
  max_year <- 4L
  n_cov <- 2L

  sim <- simulate_crirt_data(
    N_STUDENTS = n_students,
    N_EXAMS = n_exams,
    N_GRADES = n_grades,
    MAX_YEAR = max_year,
    N_COV = n_cov,
    SEED = seed,
    PARAMS = params,
    CHECK = check
  )

  if (!check) {
    sim$data_dims <- list(
      n_obs = n_students,
      n_exams = n_exams,
      n_grades = n_grades,
      n_cov = n_cov,
      yb = max_year
    )
    sim$yb <- max_year
  }

  fit <- list(
    mod = "full",
    fit = list(par = parList2Vec(sim$params)),
    data = sim
  )
  fit$data$data_dims <- sim$data_dims
  fit$data$labs$exams <- paste0("E", seq_len(n_exams))
  fit$data$yb <- sim$yb
  fit$data$data_dims$yb <- sim$yb
  fit
}

# cumulative grade probabilities
cumulatives <- function(pars, lat, j, n_grades) {
  vapply(
    seq_len(n_grades),
    function(g) stats::plogis(pars$irt[j, "slope"] * lat[, 1] - pars$irt[j, g]),
    numeric(nrow(lat))
  )
}

test_that("the grade moments match a brute force sum over the categories", {
  fit <- make_fit()
  np <- fit$data$data_dims
  lat <- fit$data$latent
  pars <- parVec2List(
    fit$fit$par,
    N_GRADES = np$n_grades,
    N_EXAMS = np$n_exams,
    N_COV = np$n_cov,
    YB = np$yb
  )

  res <- conditioned_residuals(fit, lat)

  for (j in seq_len(np$n_exams)) {
    cum <- cumulatives(pars, lat, j, np$n_grades)
    # category probabilities from the cumulatives
    p <- cum - cbind(cum[, -1, drop = FALSE], 0)
    p <- p / cum[, 1]
    x <- seq_len(np$n_grades)

    want_mean <- drop(p %*% x)
    want_var <- drop(p %*% x^2) - want_mean^2

    expect_equal(
      unname(res$grade[, j]),
      unname(fit$data$gradesMat[, j] - want_mean),
      tolerance = 1e-10
    )
    expect_equal(unname(res$grade_sd[, j]), unname(sqrt(want_var)), tolerance = 1e-10)
  }

  cum <- cumulatives(pars, lat, 1L, np$n_grades)
  p <- (cum - cbind(cum[, -1, drop = FALSE], 0)) / cum[, 1]
  expect_equal(rowSums(p), rep(1, np$n_obs), tolerance = 1e-12)
})

test_that("the time moments match numerical integration of the truncated normal", {
  fit <- make_fit()
  np <- fit$data$data_dims
  lat <- fit$data$latent
  pars <- parVec2List(
    fit$fit$par,
    N_GRADES = np$n_grades,
    N_EXAMS = np$n_exams,
    N_COV = np$n_cov,
    YB = np$yb
  )

  res <- conditioned_residuals(fit, lat)

  for (j in c(1L, 3L, np$n_exams)) {
    mu <- unname(pars$irt[j, "time_loc"] - lat[, 2])
    sd_time <- unname(1 / pars$irt[j, "time_invsd"])
    upper <- unname(log(fit$data$max_time))

    for (i in c(1L, 7L, 25L, np$n_obs)) {
      mass <- stats::pnorm(upper[i], mu[i], sd_time)
      m1 <- stats::integrate(
        function(z) z * stats::dnorm(z, mu[i], sd_time),
        -Inf,
        upper[i]
      )$value /
        mass
      # variance taken about the mean
      v <- stats::integrate(
        function(z) (z - m1)^2 * stats::dnorm(z, mu[i], sd_time),
        -Inf,
        upper[i]
      )$value /
        mass

      expect_equal(
        unname(res$time[i, j]),
        unname(log(fit$data$timeMat[i, j]) - m1),
        tolerance = 1e-6
      )
      expect_equal(unname(res$time_sd[i, j]), sqrt(v), tolerance = 1e-6)
    }
  }
})

test_that("the four matrices are shaped alike and carry the right missingness", {
  fit <- make_fit()
  np <- fit$data$data_dims
  res <- conditioned_residuals(fit, fit$data$latent)

  expect_named(res, c("grade", "time", "grade_sd", "time_sd"))
  for (m in res) {
    expect_equal(dim(m), c(np$n_obs, np$n_exams))
    expect_equal(colnames(m), fit$data$labs$exams)
  }

  # numerators missing where data missing
  expect_identical(unname(is.na(res$grade)), unname(is.na(fit$data$gradesMat)))
  expect_identical(unname(is.na(res$time)), unname(is.na(fit$data$timeMat)))

  expect_false(anyNA(res$grade_sd))
  expect_false(anyNA(res$time_sd))
  expect_true(any(is.na(res$grade) & !is.na(res$grade_sd)))
})

test_that("LAST_DAY defaults to the last observed day and drives only the times", {
  fit <- make_fit()
  res <- conditioned_residuals(fit, fit$data$latent)

  expect_identical(
    res,
    conditioned_residuals(fit, fit$data$latent, LAST_DAY = fit$data$max_time)
  )

  # later truncation lowers the residual
  later <- conditioned_residuals(
    fit,
    fit$data$latent,
    LAST_DAY = fit$data$max_time * 10
  )
  expect_identical(later$grade, res$grade)
  expect_identical(later$grade_sd, res$grade_sd)

  keep <- !is.na(res$time)
  expect_true(all(later$time[keep] < res$time[keep]))
  expect_true(all(later$time_sd > res$time_sd))

  expect_error(conditioned_residuals(fit, fit$data$latent, LAST_DAY = 1))
})

test_that("standardised residuals are centred and unit scaled when censoring is exogenous", {
  # close every exit fixing censoring
  base <- make_fit()
  np <- base$data$data_dims
  pars <- parVec2List(
    base$fit$par,
    N_GRADES = np$n_grades,
    N_EXAMS = np$n_exams,
    N_COV = np$n_cov,
    YB = np$yb
  )
  pars$cr$beta[seq_len(np$yb), ] <- -500
  pars$cr$beta[c("ability", "speed"), ] <- 0
  pars$cr$grad <- -500

  fit <- make_fit(seed = 4, n_students = 2000, params = pars, check = FALSE)
  expect_true(all(fit$data$last_year == np$yb))
  expect_true(all(fit$data$max_time == np$yb * 365))

  res <- conditioned_residuals(fit, fit$data$latent)
  grade_std <- res$grade / res$grade_sd
  time_std <- res$time / res$time_sd

  expect_lt(abs(mean(grade_std, na.rm = TRUE)), 0.05)
  expect_lt(abs(stats::sd(grade_std, na.rm = TRUE) - 1), 0.05)
  expect_lt(abs(mean(time_std, na.rm = TRUE)), 0.05)
  expect_lt(abs(stats::sd(time_std, na.rm = TRUE) - 1), 0.05)
})

test_that("simulate_residuals() stacks the four families and records its arguments", {
  fit <- make_fit()
  np <- fit$data$data_dims

  sim <- simulate_residuals(fit, B = 3, NCORES = 1L, SEED = 2)

  expect_named(sim, c("observed", "simulated", "meta"))
  expect_named(sim$simulated, c("grade", "time", "grade_sd", "time_sd"))
  for (arr in sim$simulated) {
    expect_equal(dim(arr), c(np$n_obs, np$n_exams, 3L))
    expect_equal(dimnames(arr)[[2]], fit$data$labs$exams)
  }
  expect_identical(sim$observed, conditioned_residuals(fit))
  expect_null(sim$meta$centre_time)
  expect_identical(sim$meta$B, 3)

  # truncation varies across replications
  expect_false(identical(
    sim$simulated$time_sd[, , 1L],
    sim$simulated$time_sd[, , 2L]
  ))
})
