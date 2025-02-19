test_that("internal_lat() loglikelihood", {
  skip_if_not_installed("mvtnorm")
  set.seed(123)
  n_cov <- 5
  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  L <- t(chol(S))
  covariates <- rnorm(n_cov)
  B <- matrix(rnorm(2*n_cov), 2, n_cov)
  mu <- B%*%covariates

  pars <- c(L[2,1], L[2,2], as.numeric(t(B)))


  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    ll <- cpp_lat(
      THETA = pars,
      COVARIATES = covariates,
      ABILITY = lat[1],
      SPEED = lat[2])$ll
    Rval <- mvtnorm::dmvnorm(lat, mean = mu, sigma = S, log = TRUE)
    expect_true(abs(ll-Rval)<1e-6)
  }
})

test_that("internal_lat() gradient", {
  skip_if_not_installed("numDeriv")
  set.seed(123)
  n_cov <- 5
  sig <- 2
  rho <- .6
  S <- matrix(c(1,rho*sig,rho*sig,sig^2), 2, 2)
  L <- t(chol(S))
  covariates <- rnorm(n_cov)
  B <- matrix(rnorm(2*n_cov), 2, n_cov)
  mu <- B%*%covariates

  pars <- c(L[2,1], L[2,2], as.numeric(t(B)))



  rfun <- function(PAR, LAT){
    cpp_lat(
      THETA = PAR,
      COVARIATES = covariates,
      ABILITY = LAT[1],
      SPEED = LAT[2])$ll
  }
  set.seed(123)
  for (i in 1:5) {
    lat <- rnorm(2)
    grll <- cpp_lat(
      THETA = pars,
      COVARIATES = covariates,
      ABILITY = lat[1],
      SPEED = lat[2])$gr
    Rval <- numDeriv::grad(func = rfun, x = pars, LAT = lat)
    expect_equal(grll, Rval)
  }
})


# #############################
# ### gen params ####
# n_grades <- 4L
# n_exams <- 5L
# n_cov <- 1L
# yb <- 3
# dim_cr <- 2*(yb+n_cov+2)+1
# dim_irt_lat <- n_exams*(n_grades+3)+2
# labs_exams <- paste0('ECO0',1:n_exams)
# labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
# theta <- c(rnorm(dim_irt_lat),rep(NA, dim_cr))
# parList <- parVec2List(
#   THETA = theta,
#   N_GRADES = n_grades,
#   N_EXAMS = n_exams,
#   N_COV = n_cov,
#   YB = yb,
#   LABS_EXAMS = labs_exams,
#   LABS_GRADES = labs_grades)
# Lchol <- diag(1,2,2); Lchol[2,] <- theta[(dim_irt_lat-1):(dim_irt_lat)]
# RLAT <- function(x, LAT, ROTATE, GRFLAG){
#   if(ROTATE){
#     L <- diag(1,2,2); L[2,] <- x[(dim_irt_lat-1):(dim_irt_lat)]
#     internal_lat <- as.numeric(L%*%LAT)
#   }else{
#     internal_lat <- LAT
#   }
#
#   obj <- cpp_lat(
#     THETA = x,
#     ABILITY = internal_lat[1],
#     SPEED = internal_lat[2],
#     DIM_IRT = dim_irt_lat-2,
#     ROTATED = ROTATE,
#     GRFLAG = GRFLAG
#   )
#
#   if(GRFLAG){
#     return(obj$grll[1:dim_irt_lat])
#   }else{
#     return(obj$ll)
#   }
# }
#
# lat <- c(1,.25)
# RLAT(x=theta, LAT=lat, ROTATE = TRUE, GRFLAG = FALSE)
# numDeriv::grad(x=theta[1:dim_irt_lat], func=RLAT, LAT=lat, ROTATED = FALSE, GRFLAG = FALSE)
# test_that("check gradient without rotation", {
#   expect_equal(
#     RLAT(x=theta, LAT=lat, ROTATE = FALSE, GRFLAG = TRUE),
#     numDeriv::grad(x=theta[1:dim_irt_lat], func=RLAT, LAT=lat, ROTATE = FALSE, GRFLAG = FALSE)
#   )
# })
#
# test_that("check gradient with rotation", {
#   expect_equal(
#     RLAT(x=theta, LAT=lat, ROTATE = TRUE, GRFLAG = TRUE),
#     numDeriv::grad(x=theta[1:dim_irt_lat], func=RLAT, LAT=lat, ROTATE = TRUE, GRFLAG = FALSE)
#   )
# })
#
#
#
#
#
#
#
#
#
#
#
