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
theta_lat <- rnorm(dim_lat); theta_lat[2] <- abs(theta_lat[2])
theta_cr <- rnorm(dim_cr)
theta <- c(theta_irt, theta_lat, theta_cr)

X <- rnorm(n_cov)
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  N_COV = n_cov,
  YB = yb,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)

mat <- irtVec2Mat(
  THETA_IRT = theta_irt,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades
)




#### test exam-specific likelihood ####
FUNEXAM <- function(x, OBSFLAG, LATPARFLAG=FALSE, OUT="ll"){
  ab <- ability
  sp <- speed
  if(LATPARFLAG){
    ab <- ab + x[(dim_irt+3):(dim_irt+2+n_cov)]%*%X
    sp <- x[n_exams*(n_grades+3)+1]*ab + x[n_exams*(n_grades+3)+2]*sp + x[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)]%*%X
  }
  obj <- cpp_examLik(
    EXAM = exam-1,
    GRADE = grade,
    DAY = day,
    MAX_DAY = day,
    OBSFLAG = OBSFLAG,
    THETA = x,
    COVARIATES = X,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    ABILITY = ab,
    SPEED = sp,
    LATPARFLAG=LATPARFLAG
  )
  if(OUT=="ll"){
    return(obj$ll)
  }else{
    return(obj$grll)
  }
}



set.seed(333)
for (ability in rnorm(3,0,1)) {
  for (speed in rnorm(3,0,1)) {
    for (day in sample(100:1000, 3, replace = TRUE)) {
      for (exam in 1:n_exams) {
        for (grade in 1:n_grades) {
          pG <- cpp_pGrade(GRADE = grade, EXAM = exam-1, THETA = theta, COVARIATES = X, N_GRADES = n_grades, N_EXAMS = n_exams, ABILITY = ability, LOGFLAG = FALSE, LATPARFLAG = FALSE)$prob
          pT <- dlnorm(day,
                       mat[exam,n_grades+2]-speed,
                       1/mat[exam,n_grades+3])
          Rval <- pT*pG
          val <- FUNEXAM(
            OBSFLAG = T,
            x = theta
          )
          test_that("examLik() log output observed exam", {
            skip_if(!is.finite(log(Rval)))
            expect_equal(val, log(Rval))
          })

          pG <- cpp_pGreaterGrades(GRADE = 1, EXAM = exam-1, THETA = theta, COVARIATES = X, N_GRADES = n_grades, N_EXAMS = n_exams, ABILITY = ability, LOGFLAG = FALSE, LATPARFLAG = FALSE)$prob
          pT <- plnorm(day,
                       mat[exam, n_grades+2]-speed,
                       1/mat[exam, n_grades+3])
          Rval <- 1-pT*pG
          val <- FUNEXAM(
            OBSFLAG = F,
            x = theta
          )
          test_that("examLik() log output not observed exam", {
            skip_if(!is.finite(log(Rval)))
            expect_equal(val, log(Rval))
          })

        }

      }
    }
  }
}


## gradient checks ######
n <- 10
latMat <- matrix(c(rnorm((n-1)*2), -5,5),n,2, byrow = TRUE)
X <- matrix(rnorm(n*n_cov), n, n_cov)
set.seed(123)

#### sim grades ####
gradesMat <- matrix(0, n, n_exams)

for (i in 1:n) {
  for (e in 1:n_exams) {
    #linear predictor exams X grades
    linp <- mat[e, n_grades+1] * latMat[i,1] - mat[e, 1:n_grades]

    # probabilities of greater grades
    pgg <- exp(linp)/(1+exp(linp))

    # probabilities of grades
    pg <- c(pgg[1:(n_grades-1)] - pgg[2:n_grades], pgg[n_grades])
    gradesMat[i,e] <- which(rmultinom(n = 1, size=1, prob = c(1-sum(pg), pg))==1)-1

  }
}



#### sim times ####
set.seed(123)
timeMat <- matrix(NA, n, n_exams)
for (i in 1:n) {
  for (e in 1:n_exams) {
    timeMat[i,e] <- exp(
      rnorm(1,
            mean = mat[e, n_grades + 2]-latMat[i,2],
            sd = 1/mat[e, n_grades + 3])+1
    )
  }
}
timeMat[gradesMat==0] <- NA
timeMat <- ceiling(timeMat)
#### mat to-do ####
todoMat <- matrix(1, n, n_exams)

#### censoring ####
max_day <- as.integer(max(timeMat, na.rm = TRUE)+10)
timeMat[timeMat>max_day] <- NA
obsMat <- matrix(1, n, n_exams)
obsMat[is.na(timeMat)] <- 0

#### checks ####
RFUN <- function(x, SPEED, ABILITY, LATPARFLAG, COVARIATES, OUT="ll",  ...){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(x[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- x[dim_irt+1]*ABILITY + x[dim_irt+2]*SPEED + t(x[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
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

for (i in 1:n) {
  for (exam in 1:n_exams) {
    test_that(paste0("exam density not rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")

      numGrad <- numDeriv::grad(func = RFUN, x = theta[1:(dim_irt+dim_lat)], MAX_DAY = max_day,
                                EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
                                SPEED = latMat[i,2], ABILITY = latMat[i,1], LATPARFLAG = FALSE,
                                OBSFLAG = obsMat[i, exam], COVARIATES= X[i,])
      grcpp <- RFUN(
        x = theta, MAX_DAY = max_day,
        EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
        SPEED = latMat[i,2], ABILITY = latMat[i,1], LATPARFLAG = FALSE,
        OBSFLAG = obsMat[i, exam], COVARIATES= X[i,],
        OUT="grll"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })

    test_that(paste0("exam density rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")

      numGrad <- numDeriv::grad(func = RFUN, x = theta[1:(dim_irt+dim_lat)], MAX_DAY = max_day,
                                EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
                                SPEED = latMat[i,2], ABILITY = latMat[i,1], LATPARFLAG = TRUE,
                                OBSFLAG = obsMat[i, exam], COVARIATES= X[i,])
      grcpp <- RFUN(
        x = theta, MAX_DAY = max_day,
        EXAM = exam-1, DAY = timeMat[i, exam], GRADE = gradesMat[i, exam],
        SPEED = latMat[i,2], ABILITY = latMat[i,1], LATPARFLAG = TRUE,
        OBSFLAG = obsMat[i, exam], COVARIATES= X[i,],
        OUT="grll"
      )

      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })


  }

}

