n_grades <- 4L
n_exams  <- 3L
n_cov    <- 10L
yb       <- 5L

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



#### test times probabilities output values #####
set.seed(111)
speeds <- sort(rnorm(3,0,2), decreasing=T)
FUNTIME <- function(x, CDFFLAG, LOGFLAG, SPEED, ABILITY=0, ROTATED=FALSE, OUT="prob"){
  obj <- cpp_pTimeExam(
    EXAM = exam-1,
    DAY = day,
    THETA = x,
    COVARIATES = X,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    SPEED = SPEED,
    ABILITY = ABILITY,
    CDFFLAG = CDFFLAG,
    LOGFLAG = LOGFLAG,
    LATPARFLAG = ROTATED
  )
  if(OUT=="prob"){
    return(obj$prob)
  }else{
    return(obj$gr[1:dim_irt_lat])
  }
}
for (speed_index in 1:length(speeds)) {
  for (day in runif(10, 100, 1000)) {
    for (exam in 1:n_exams) {



      Rval <- dlnorm(day,
                     mat[exam, n_grades+2]-speeds[speed_index],
                     1/mat[exam, n_grades+3])
      val <- FUNTIME(
        x = theta,
        CDFFLAG = FALSE,
        LOGFLAG = FALSE,
        SPEED = speeds[speed_index]
      )

      test_that(paste0("pTimeExam() val, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam), {
        expect_equal(val, Rval)
      })

      Rval <- plnorm(day,
                     mat[exam, n_grades+2]-speeds[speed_index],
                     1/mat[exam, n_grades+3])
      val <- FUNTIME(
        x = theta,
        CDFFLAG = TRUE,
        LOGFLAG = FALSE,
        SPEED = speeds[speed_index]
      )

      test_that(paste0("pTimeExam() cdf, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  expect_equal(val, Rval)
                })


      if(speed_index>1){
        val <- FUNTIME(
          x = theta,
          CDFFLAG = TRUE,
          LOGFLAG = FALSE,
          SPEED = speeds[speed_index]
        )
        valprev <- FUNTIME(
          x = theta,
          CDFFLAG = TRUE,
          LOGFLAG = FALSE,
          SPEED = speeds[speed_index-1]
        )
        # test_that(paste0("pTimeExam() val decreases with decresing speeds, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
        #           {
        #             expect_true(valprev>val)
        #           })
      }




    }
  }

}

speeds <- sort(c(1000, -1000), decreasing=T)

set.seed(222)
for (speed_index in 1:length(speeds)) {
  for (day in runif(3, 100, 1000)) {
    for (exam in 1:n_exams) {


      val <- FUNTIME(
        x = theta,
        SPEED = speeds[speed_index],
        CDFFLAG = F,
        LOGFLAG = TRUE
      )

      test_that(paste0("pTimeExam() log extreme, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  expect_true(is.finite(val))
                })

      val <- FUNTIME(
        x = theta,
        SPEED = speeds[speed_index],
        CDFFLAG = T,
        LOGFLAG = TRUE
      )
      test_that(paste0("pTimeExam() cdf log extreme, speed:", round(speeds[speed_index],2), ", day:", round(day,2), ", exam:", exam),
                {
                  expect_true(is.finite(val))
                })
    }
  }

}

########## gradient tests ##########

n <- 10
n_grades <- 4L
n_exams  <- 3L
n_cov    <- 10L
yb       <- 5L

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
latMat <- matrix(c(rnorm((n-1)*2), -5,5),n,2, byrow = TRUE)
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
timeMat <- ceiling(timeMat)


#### checks #####

FUNTIME <- function(x, EXAM, DAY, CDFFLAG, LOGFLAG, SPEED, COVARIATES=X, ABILITY=0, LATPARFLAG=FALSE, OUT="prob"){
  ab <- ABILITY
  sp <- SPEED
  if(LATPARFLAG){
    ab <- ABILITY + t(x[(dim_irt+3):(dim_irt+2+n_cov)])%*%COVARIATES
    sp <- x[dim_irt+1]*ABILITY + exp(x[dim_irt+2])*SPEED + t(x[(dim_irt+2+n_cov+1):(dim_irt+dim_lat)])%*%COVARIATES
  }

  obj <- cpp_pTimeExam(
    EXAM = EXAM,
    DAY = DAY,
    THETA=x,
    COVARIATES = COVARIATES,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    SPEED = sp,
    ABILITY = ab,
    CDFFLAG = CDFFLAG,
    LOGFLAG = LOGFLAG,
    LATPARFLAG = LATPARFLAG
  )
  if(OUT=="prob"){
    return(obj$prob)
  }else{
    return(obj$gr[1:(dim_irt+dim_lat)])
  }
}

for (i in 1:n) {
  for (exam in 1:n_exams) {
    test_that(paste0("gradient time density non rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta[1:(dim_irt+dim_lat)],
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = FALSE, LATPARFLAG = FALSE)
      grcpp <- FUNTIME(
        x = theta,
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        SPEED = latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        LATPARFLAG = FALSE,
        LOGFLAG = FALSE,
        OUT="gr"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })
    test_that(paste0("gradient time log density non rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta[1:(dim_irt+dim_lat)],
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = TRUE, LATPARFLAG = FALSE)
      grcpp <- FUNTIME(
        x = theta,
        EXAM = exam-1, DAY = timeMat[i, exam],
        SPEED = latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        LATPARFLAG = FALSE,
        LOGFLAG = TRUE,
        OUT="gr"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })

    test_that(paste0("gradient time density rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta[1:(dim_irt+dim_lat)],
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = FALSE, LATPARFLAG = TRUE)
      grcpp <- FUNTIME(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        x = theta,
        SPEED = latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        LATPARFLAG = TRUE,
        LOGFLAG = FALSE,
        OUT="gr"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)

    })
    test_that(paste0("gradient time log density rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta[1:(dim_irt+dim_lat)],
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = FALSE, LOGFLAG = TRUE,
                                LATPARFLAG = TRUE)
      grcpp <- FUNTIME(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        x = theta,
        SPEED = latMat[i,2],
        CDFFLAG = FALSE,
        ABILITY = latMat[i,1],
        LATPARFLAG = TRUE,
        LOGFLAG = TRUE,
        OUT="gr"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)

    })

    test_that(paste0("gradient time cdf non rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta[1:(dim_irt+dim_lat)],
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = TRUE, LOGFLAG = FALSE, LATPARFLAG = FALSE)
      grcpp <- FUNTIME(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        x = theta,
        SPEED = latMat[i,2],
        CDFFLAG = TRUE,
        ABILITY = latMat[i,1],
        LATPARFLAG = FALSE,
        LOGFLAG = FALSE,
        OUT="gr"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)
    })

    test_that(paste0("gradient time cdf rotated lat, i:",i, ", e:",exam), {
      skip_if_not_installed("numDeriv")
      skip_if(is.na(timeMat[i, exam]))

      numGrad <- numDeriv::grad(func = FUNTIME, x = theta[1:(dim_irt+dim_lat)],
                                EXAM = exam-1, DAY = timeMat[i, exam], SPEED = latMat[i,2],
                                ABILITY = latMat[i,1], CDFFLAG = TRUE, LOGFLAG = FALSE, LATPARFLAG = TRUE)
      grcpp <- FUNTIME(
        EXAM = exam-1,
        DAY = timeMat[i, exam],
        x = theta,
        SPEED = latMat[i,2],
        CDFFLAG = TRUE,
        ABILITY = latMat[i,1],
        LATPARFLAG = TRUE,
        LOGFLAG = FALSE,
        OUT="gr"
      )
      expect_equal(grcpp, numGrad,  tolerance = 1e-5)

    })
  }
}





