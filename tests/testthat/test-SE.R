# pin the delta method orientation
n_grades <- 4L
n_exams <- 3L
n_cov <- 3L
yb <- 5L

dim_irt <- n_exams * (n_grades + 3)
dim_lat <- 2 + 2 * n_cov
dim_cr <- 2 * (yb + 2) + 1

labs_exams <- paste0("ECO0", 1:n_exams)
labs_grades <- c("[18,22)", "[22,25)", "[25,28)", "[29,30L]")
labs_cov <- paste0("X", 1:n_cov)

set.seed(123)
theta_lat <- rnorm(dim_lat)
theta_lat[2] <- log(abs(theta_lat[2]))
theta <- c(rnorm(dim_irt), theta_lat, rnorm(dim_cr))

repar <- function(x) {
  parVec2Repar(
    THETA = x,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    N_COV = n_cov,
    YB = yb,
    LABS_EXAMS = labs_exams,
    LABS_GRADES = labs_grades,
    LABS_COV = labs_cov
  )
}

# small covariance keeps linearisation accurate
p <- length(theta)
V <- crossprod(matrix(rnorm(p * p), p, p)) / p * 1e-4

jacob <- numDeriv::jacobian(func = repar, x = theta)
se_delta <- sqrt(rowSums((jacob %*% V) * jacob))
se_transposed <- sqrt(diag(t(jacob) %*% V %*% jacob))

# reference distribution by simulation
draws <- sweep(matrix(rnorm(2e4 * p), ncol = p) %*% chol(V), 2, theta, "+")
se_mc <- apply(apply(draws, 1, repar), 1, stats::sd)

tidy <- parVec2Repar(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  N_COV = n_cov,
  YB = yb,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades,
  LABS_COV = labs_cov,
  TIDY = TRUE
)

max_rel <- function(x, y) max(abs(x / y - 1))

test_that("the delta method reproduces the sampling distribution", {
  expect_lt(max_rel(se_delta, se_mc), 0.04)
})

test_that("the transposed product does not", {
  expect_gt(max_rel(se_transposed, se_mc), 0.5)
})

test_that("orientation bites on the thresholds and the correlation", {
  # nonlinear entries expose the orientation
  thresholds <- which(tidy$type %in% labs_grades)
  correlation <- which(tidy$type == "correlation")
  expect_length(correlation, 1)
  expect_gt(length(thresholds), 0)
  expect_gt(max_rel(se_transposed[thresholds], se_delta[thresholds]), 0.5)
  expect_gt(max_rel(se_transposed[correlation], se_delta[correlation]), 0.2)
})

test_that("parameters copied through are unaffected by the orientation", {
  # copied entries hide the orientation
  copied <- which(tidy$type %in% c("slope", "time_loc"))
  expect_gt(length(copied), 0)
  expect_lt(max_rel(se_transposed[copied], se_delta[copied]), 1e-8)
})
