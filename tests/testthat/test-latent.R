test_that("internal_lat() loglikelihood", {
  skip_if_not_installed("mvtnorm")
  set.seed(123)
  n_cov <- 5
  sig <- 2
  rho <- .6
  S <- matrix(c(1, rho * sig, rho * sig, sig^2), 2, 2)
  L <- t(chol(S))
  covariates <- rnorm(n_cov)
  B <- matrix(rnorm(2 * n_cov), 2, n_cov)
  mu <- B %*% covariates

  pars <- c(L[2, 1], log(L[2, 2]), as.numeric(t(B)))

  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    ll <- cpp_lat(
      THETA = pars,
      COVARIATES = covariates,
      ABILITY = lat[1],
      SPEED = lat[2]
    )$ll
    Rval <- mvtnorm::dmvnorm(lat, mean = mu, sigma = S, log = TRUE)
    expect_true(abs(ll - Rval) < 1e-6)
  }
})

test_that("internal_lat() gradient", {
  skip_if_not_installed("numDeriv")
  set.seed(123)
  n_cov <- 5
  sig <- 2
  rho <- .6
  S <- matrix(c(1, rho * sig, rho * sig, sig^2), 2, 2)
  L <- t(chol(S))
  covariates <- rnorm(n_cov)
  B <- matrix(rnorm(2 * n_cov), 2, n_cov)
  mu <- B %*% covariates

  pars <- c(L[2, 1], log(L[2, 2]), as.numeric(t(B)))

  rfun <- function(PAR, LAT) {
    cpp_lat(
      THETA = PAR,
      COVARIATES = covariates,
      ABILITY = LAT[1],
      SPEED = LAT[2]
    )$ll
  }
  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    grll <- cpp_lat(
      THETA = pars,
      COVARIATES = covariates,
      ABILITY = lat[1],
      SPEED = lat[2]
    )$grll
    Rval <- numDeriv::grad(func = rfun, x = pars, LAT = lat)
    expect_equal(grll, Rval)
  }
})
