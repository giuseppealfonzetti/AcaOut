set.seed(123)

### gen params ####
n_grades <- 2L
n_exams  <- 2L
n_cov    <- 1
yb       <- 3L

dim_irt <- n_exams*(n_grades+3)
dim_lat <- 2+2*n_cov
dim_cr  <- 2*(yb+2)+1

labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)')
labs_cov <- if(n_cov>0) paste0("X",1:n_cov)

set.seed(123)
theta_irt <- rnorm(dim_irt)
theta_lat <- rnorm(dim_lat); theta_lat[2] <- log(abs(theta_lat[2]))
theta_cr <- rnorm(dim_cr)
theta <- c(theta_irt, theta_lat, theta_cr)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  N_COV = n_cov,
  YB = yb,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)


### hazard checks
RFUN <- function(PAR, ABILITY, SPEED, YEAR, OUTCOME, COVARIATES, LATPARFLAG, GRFLAG=FALSE){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(PAR[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- PAR[dim_irt+1]*ABILITY + exp(PAR[dim_irt+2])*SPEED + t(PAR[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }
  obj <- cpp_hazard(
    OUTCOME=OUTCOME,
    YEAR = YEAR,
    THETA = PAR,
    COVARIATES = COVARIATES,
    ABILITY = ab,
    SPEED = sp,
    YB = yb,
    LATPARFLAG = LATPARFLAG)

  if(GRFLAG){return(obj$gr)}else{obj$prob}
}

lat <- rnorm(2)
covariates <- rnorm(n_cov)
for (year in 1:yb) {
  for (outcome in 1:3) {
    test_that(paste0("Hazard year ", year, ", outcome ",outcome),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=FALSE, GRFLAG=TRUE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=FALSE, GRFLAG=FALSE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)]
      )
    })

    test_that(paste0("Hazard year ", year, ", outcome ",outcome, " | LATPARFLAG=TRUE"),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=TRUE, GRFLAG=TRUE)[(1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=TRUE, GRFLAG=FALSE)[(1):(dim_irt+dim_lat+dim_cr)]
      )
    })


  }
}

### survival checks
RFUN <- function(PAR, ABILITY, SPEED, YEAR, COVARIATES, LATPARFLAG, GRFLAG=FALSE){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(PAR[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- PAR[dim_irt+1]*ABILITY + exp(PAR[dim_irt+2])*SPEED + t(PAR[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }
  obj <- cpp_survival(
    YEAR_FIRST = 1,
    YEAR_LAST = YEAR,
    THETA = PAR,
    COVARIATES = COVARIATES,
    ABILITY = ab,
    SPEED = sp,
    YB = yb,
    YEAR_LAST_EXAM = 10,
    LATPARFLAG = LATPARFLAG)

  if(GRFLAG){return(obj$gr)}else{obj$prob}
}


lat <- rnorm(2)
covariates <- rnorm(n_cov)
for (year in 1:yb) {
    test_that(paste0("Survival year ", year, ", outcome ",outcome),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, COVARIATES=covariates, LATPARFLAG=FALSE, GRFLAG=TRUE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year,  COVARIATES=covariates, LATPARFLAG=FALSE, GRFLAG=FALSE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)]
      )
    })

    test_that(paste0("Survival year ", year, ", outcome ",outcome, " | LATPARFLAG=TRUE"),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year,  COVARIATES=covariates, LATPARFLAG=TRUE, GRFLAG=TRUE)[(1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year,  COVARIATES=covariates, LATPARFLAG=TRUE, GRFLAG=FALSE)[(1):(dim_irt+dim_lat+dim_cr)]
      )
    })
}

### outcome checks
RFUN <- function(PAR, ABILITY, SPEED, YEAR, OUTCOME, COVARIATES, YLE, LATPARFLAG, GRFLAG=FALSE){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(PAR[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- PAR[dim_irt+1]*ABILITY + exp(PAR[dim_irt+2])*SPEED + t(PAR[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }
  obj <- cpp_outcome(
    OUTCOME = OUTCOME,
    YEAR_FIRST = 1,
    YEAR_LAST = YEAR,
    THETA = PAR,
    COVARIATES = COVARIATES,
    ABILITY = ab,
    SPEED = sp,
    YB = yb,
    YEAR_LAST_EXAM = YLE,
    LATPARFLAG = LATPARFLAG)

  if(GRFLAG){return(obj$grl)}else{obj$logprob}
}

lat <- rnorm(2)
covariates <- rnorm(n_cov)
for (year in 1:yb) {
  for (outcome in 0:3) {
    test_that(paste0("Survival year ", year, ", outcome ",outcome),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=FALSE, YLE=year+1, GRFLAG=TRUE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome,  COVARIATES=covariates, YLE=year+1, LATPARFLAG=FALSE, GRFLAG=FALSE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)]
      )
    })

    test_that(paste0("Survival year ", year, ", outcome ",outcome, " | LATPARFLAG=TRUE"),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=TRUE, YLE=year+1, GRFLAG=TRUE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome,  COVARIATES=covariates, YLE=year+1, LATPARFLAG=TRUE, GRFLAG=FALSE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)]
      )
    })
  }
}

### ccr class checks
RFUN <- function(PAR, ABILITY, SPEED, YEAR, OUTCOME, COVARIATES, YLE, LATPARFLAG, GRFLAG=FALSE){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(PAR[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- PAR[dim_irt+1]*ABILITY + exp(PAR[dim_irt+2])*SPEED + t(PAR[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }
  obj <- cpp_cr_class(
    OUTCOME = OUTCOME,
    YEAR_FIRST = 1,
    YEAR_LAST = YEAR,
    THETA = PAR,
    COVARIATES = COVARIATES,
    ABILITY = ab,
    SPEED = sp,
    YB = yb,
    YEAR_LAST_EXAM = YLE,
    LATPARFLAG = LATPARFLAG)

  if(GRFLAG){return(obj$grll)}else{obj$ll}
}

lat <- rnorm(2)
covariates <- rnorm(n_cov)
for (year in 1:yb) {
  for (outcome in 0:3) {
    test_that(paste0("Survival year ", year, ", outcome ",outcome),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=FALSE, YLE=year+1, GRFLAG=TRUE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome,  COVARIATES=covariates, YLE=year+1, LATPARFLAG=FALSE, GRFLAG=FALSE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)]
      )
    })

    test_that(paste0("Survival year ", year, ", outcome ",outcome, " | LATPARFLAG=TRUE"),{
      skip_if_not_installed("numDeriv")
      expect_equal(
        RFUN(PAR = theta, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome, COVARIATES=covariates, LATPARFLAG=TRUE, YLE=year+1, GRFLAG=TRUE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)],
        numDeriv::grad(x=theta, func = RFUN, ABILITY=lat[1], SPEED=lat[2], YEAR=year, OUTCOME=outcome,  COVARIATES=covariates, YLE=year+1, LATPARFLAG=TRUE, GRFLAG=FALSE)[(dim_irt+1):(dim_irt+dim_lat+dim_cr)]
      )
    })
  }
}

