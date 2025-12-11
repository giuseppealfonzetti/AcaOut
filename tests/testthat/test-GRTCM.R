n <- 10
set.seed(123)

### gen params ####
n_grades <- 4L
n_exams  <- 3L
n_cov    <- 2
yb       <- 2L

dim_irt <- n_exams*(n_grades+3)
dim_lat <- 2+2*n_cov
dim_cr  <- 2*(yb+2)+1

labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
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
irtMat <- parList$irt
# theta_irt <- irtMat2Vec(irtMat)
X <- matrix(rnorm(n*n_cov), n, n_cov)
latMat <- matrix(rnorm(n*2), n, 2)
nodes <- expand.grid(c(1,-1),c(1,-1))
weights <- rep(.25, 4)

#### sim grades ####
gradesMat <- matrix(0, n, n_exams)

for (i in 1:n) {
  for (e in 1:n_exams) {
    #linear predictor exams X grades
    linp <- irtMat[e, n_grades+1] * latMat[i,1] - irtMat[e, 1:n_grades]

    # probabilities of greater grades
    pgg <- exp(linp)/(1+exp(linp))

    # probabilities of grades
    pg <- c(pgg[1:(n_grades-1)] - pgg[2:n_grades], pgg[n_grades])
    gradesMat[i,e] <- which(rmultinom(n = 1, size=1, prob = c(1-sum(pg), pg))==1)-1

  }
}



#### sim times ####
set.seed(123)
timeMat <- matrix(0, n, n_exams)
for (i in 1:n) {
  for (e in 1:n_exams) {
    timeMat[i,e] <- exp(
      rnorm(1,
            mean = irtMat[e, n_grades + 2]-latMat[i,2],
            sd = 1/irtMat[e, n_grades + 3])
    )
  }
}
timeMat[gradesMat==0] <- NA

#### mat to-do ####
todoMat <- matrix(1, n, n_exams)

#### censoring ####
max_day <- max(timeMat, na.rm = TRUE)+10
timeMat[timeMat>max_day] <- NA
obsMat <- matrix(1, n, n_exams)
obsMat[is.na(timeMat)] <- 0


##### TEST GRTCM CLASS #####
RFUN <- function(x, ID, COVARIATES=X[ID,], LATPARFLAG, LAT_POINTS, GRFLAG=FALSE){
  ab <- LAT_POINTS[1]
  sp <- LAT_POINTS[2]
  if(LATPARFLAG){
    ab <- LAT_POINTS[1] + t(x[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- x[dim_irt+1]*LAT_POINTS[1] + exp(x[dim_irt+2])*LAT_POINTS[2] + t(x[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }
  obj <- cpp_grtcm_class(
    THETA = x,
    EXAMS_GRADES = gradesMat[ID,],
    EXAMS_DAYS = timeMat[ID,],
    EXAMS_SET = todoMat[ID,],
    EXAMS_OBSFLAG = obsMat[ID,],
    COVARIATES = COVARIATES,
    ABILITY = ab,
    SPEED = sp,
    MAX_DAY = max_day,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    LATPARFLAG = LATPARFLAG
  )

  if(GRFLAG){
    obj$grll
  }else{
    obj$ll
  }
}


examLik <- function(x, SPEED, ABILITY, LATPARFLAG, COVARIATES, OUT="ll",  ...){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(x[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- x[dim_irt+1]*ABILITY + exp(x[dim_irt+2])*SPEED + t(x[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }


  obj <- cpp_examLik(
    THETA = x,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    ABILITY = ab,
    SPEED = sp,
    COVARIATES = COVARIATES,
    LATPARFLAG=LATPARFLAG,

    ...
  )
  if(OUT=="ll"){
    return(obj$ll)
  }else{
    return(obj$grll[1:(dim_irt+dim_lat)])
  }
}
#### TEST GRTCM QUADRATURE ####
for (i in 1:n) {
  ll=0
  for (exam in 1:n_exams) {

    ell <- examLik(
      x = theta, MAX_DAY = max_day,
      EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
      SPEED = latMat[i,2], ABILITY = latMat[i,1], LATPARFLAG = FALSE,
      OBSFLAG = obsMat[i, exam], COVARIATES= X[i,],
      OUT="ll"
    )

    # cat("\nR| Exam", exam, ":", ell, "\n")
    ll <- ll + ell



  }


  test_that(paste0("Check student ",i, " full exam list likelihood"), {
    expect_equal(
      RFUN(x=theta, ID=i, LATPARFLAG = FALSE, LAT_POINTS = latMat[i,]),
      ll
    )



  })

  ### LATPARFLAG = TRUE
  ll=0
  for (exam in 1:n_exams) {

    ell <- examLik(
      x = theta, MAX_DAY = max_day,
      EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
      SPEED = latMat[i,2], ABILITY = latMat[i,1], LATPARFLAG = TRUE,
      OBSFLAG = obsMat[i, exam], COVARIATES= X[i,],
      OUT="ll"
    )

    # cat("\nR| Exam", exam, ":", ell, "\n")
    ll <- ll + ell



  }


  test_that(paste0("Check student ",i, " full exam list likelihood"), {
    expect_equal(
      RFUN(x=theta, ID=i, LATPARFLAG = TRUE, LAT_POINTS = latMat[i,]),
      ll
    )



  })

}

#### Check gr ####

for (i in 1:n) {
  test_that(paste0("Check student ",i, " gradient LATPARFLAG=FALSE"), {
    expect_equal(
      numDeriv::grad(RFUN, x=theta, ID=i, LATPARFLAG = FALSE, LAT_POINTS = latMat[i,]),
      RFUN(x=theta, ID=i, LATPARFLAG = FALSE, LAT_POINTS = latMat[i,], GRFLAG=TRUE)
    )
  })

  test_that(paste0("Check student ",i, " gradient LATPARFLAG=FALSE"), {
    expect_equal(
      numDeriv::grad(RFUN, x=theta, ID=i, LATPARFLAG = TRUE, LAT_POINTS = latMat[i,]),
      RFUN(x=theta, ID=i, LATPARFLAG = TRUE, LAT_POINTS = latMat[i,], GRFLAG=TRUE)
    )
  })

}

# ##### TESTS ####
# RFUN <- function(x, ROTATE){
#   GRTCM_GH(
#     THETA = x,
#     EXAMS_GRADES = gradesMat,
#     EXAMS_DAYS = timeMat,
#     EXAMS_SET = todoMat,
#     EXAMS_OBSFLAG = obsMat,
#     COVARIATES = X,
#     MAX_DAY = rep(max_day,n),
#     GRID = as.matrix(nodes),
#     WEIGHTS = weights,
#     N_GRADES = n_grades,
#     N_EXAMS = n_exams,
#     GRFLAG = FALSE,
#     ROTGRID = ROTATE
#   )$ll
# }
#
# RFUN(theta[1:(dim_irt+dim_lat)], TRUE)
#
# #### check gradient
# test_that("Check gradient derivative without node rotation", {
#
#
#   fit <- GRTCM_GH(
#     THETA = theta,
#     EXAMS_GRADES = gradesMat,
#     EXAMS_DAYS = timeMat,
#     EXAMS_SET = todoMat,
#     EXAMS_OBSFLAG = obsMat,
#     COVARIATES = X,
#     MAX_DAY = rep(max_day,n),
#     GRID = as.matrix(nodes),
#     WEIGHTS = weights,
#     N_GRADES = n_grades,
#     N_EXAMS = n_exams,
#     GRFLAG = TRUE,
#     ROTGRID = FALSE
#   )
#
#   expect_equal(
#     numDeriv::grad(RFUN, x = theta, ROTATE = FALSE),
#     fit$gr
#   )
# })
#
# test_that("Check gradient derivative with node rotation", {
#
#   fit <- GRTCM_GH(
#     THETA = theta,
#     EXAMS_GRADES = gradesMat,
#     EXAMS_DAYS = timeMat,
#     EXAMS_SET = todoMat,
#     EXAMS_OBSFLAG = obsMat,
#     COVARIATES = X,
#     MAX_DAY = rep(max_day,n),
#     GRID = as.matrix(nodes),
#     WEIGHTS = weights,
#     N_GRADES = n_grades,
#     N_EXAMS = n_exams,
#     GRFLAG = TRUE,
#     ROTGRID = TRUE
#   )
#
#   expect_equal(
#     numDeriv::grad(RFUN, x = theta, ROTATE = TRUE),
#     fit$gr,
#     tolerance = 1e-4
#   )
# })
#
##### TESTS ####
RFUN <- function(x, LATPARFLAG){
  cpp_GQ(
    THETA = x,
    EXAMS_GRADES = gradesMat,
    EXAMS_DAYS = timeMat,
    EXAMS_SET = todoMat,
    EXAMS_OBSFLAG = obsMat,
    COVARIATES = X,
    MAX_DAY = rep(max_day,n),
    OUTCOME = rep(1,n),
    YEAR_FIRST = rep(1,n),
    YEAR_LAST = rep(1,n),
    YEAR_LAST_EXAM = rep(1,n),
    YB = yb,
    GRID = as.matrix(nodes),
    WEIGHTS = weights,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    GRFLAG = FALSE,
    LATPARFLAG = LATPARFLAG,
    MOD="grtcm"
  )$ll
}

RFUN(theta, TRUE)

#### check gradient
test_that("Check gradient derivative without node rotation", {


  fit <- cpp_GQ(
    THETA = theta,
    EXAMS_GRADES = gradesMat,
    EXAMS_DAYS = timeMat,
    EXAMS_SET = todoMat,
    EXAMS_OBSFLAG = obsMat,
    COVARIATES = X,
    MAX_DAY = rep(max_day,n),
    OUTCOME = rep(1,n),
    YEAR_FIRST = rep(1,n),
    YEAR_LAST = rep(1,n),
    YEAR_LAST_EXAM = rep(1,n),
    YB = yb,
    GRID = as.matrix(nodes),
    WEIGHTS = weights,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    GRFLAG = TRUE,
    LATPARFLAG = FALSE,
    MOD="grtcm"
  )

  expect_equal(
    numDeriv::grad(RFUN, x = theta, LATPARFLAG = FALSE),
    fit$gr
  )
})

test_that("Check gradient derivative with node rotation", {


  fit <- cpp_GQ(
    THETA = theta,
    EXAMS_GRADES = gradesMat,
    EXAMS_DAYS = timeMat,
    EXAMS_SET = todoMat,
    EXAMS_OBSFLAG = obsMat,
    COVARIATES = X,
    MAX_DAY = rep(max_day,n),
    OUTCOME = rep(1,n),
    YEAR_FIRST = rep(1,n),
    YEAR_LAST = rep(1,n),
    YEAR_LAST_EXAM = rep(1,n),
    YB = yb,
    GRID = as.matrix(nodes),
    WEIGHTS = weights,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    GRFLAG = TRUE,
    LATPARFLAG = TRUE,
    MOD="grtcm"
  )

  expect_equal(
    numDeriv::grad(RFUN, x = theta, LATPARFLAG = TRUE),
    fit$gr
  )
})



######## complete loglik tests ######
RFUN <- function(x, ID, COVARIATES=X[ID,], LAT_POINTS, GRFLAG=FALSE){
  obj <- cpp_grtcm_class(
    THETA = x,
    EXAMS_GRADES = gradesMat[ID,],
    EXAMS_DAYS = timeMat[ID,],
    EXAMS_SET = todoMat[ID,],
    EXAMS_OBSFLAG = obsMat[ID,],
    COVARIATES = COVARIATES,
    ABILITY = LAT_POINTS[1],
    SPEED = LAT_POINTS[2],
    MAX_DAY = max_day,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    LATPARFLAG = FALSE
  )

  if(GRFLAG){
    obj$grcll
  }else{
    obj$cll
  }
}

for (i in 1:n) {
  test_that(paste0("Check student ",i, " gradient LATPARFLAG=FALSE"), {
    expect_equal(
      numDeriv::grad(RFUN, x=theta, ID=i, LAT_POINTS = latMat[i,]),
      RFUN(x=theta, ID=i, LAT_POINTS = latMat[i,], GRFLAG=TRUE),
      tolerance = 1e-5
    )
  })

}


