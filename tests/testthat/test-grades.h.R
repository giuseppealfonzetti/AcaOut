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




############################
RpGreaterGrades <- function(GRADE, EXAM,
                             THETA,
                             N_GRADES,
                             N_EXAMS,
                             ABILITY,
                            LOGFLAG,
                            OUT="prob",
                            LATPARFLAG){

  ab <- ABILITY
  if(LATPARFLAG & n_cov>0) ab <- ab + t(THETA[(dim_irt+3):(dim_irt+2+n_cov)])%*%X

  obj <- cpp_pGreaterGrades(
    GRADE=GRADE,
    EXAM=EXAM,
    THETA_IRT=THETA[1:dim_irt],
    THETA_LAT=THETA[(dim_irt+1):(dim_irt+dim_lat)],
    COVARIATES = X,
    N_GRADES=N_GRADES,
    N_EXAMS=N_EXAMS,
    ABILITY=ab,
    LOGFLAG=LOGFLAG,
    LATPARFLAG=LATPARFLAG)

  if(OUT=="prob"){
    return(obj$prob)
  }else{
    return(obj$gr[1:(dim_irt+dim_lat)])
  }
}

RpGrade <- function(GRADE, EXAM,
                    THETA,
                    N_GRADES,
                    N_EXAMS,
                    ABILITY,
                    LOGFLAG,
                    OUT="prob",
                    LATPARFLAG){
  ab <- ABILITY
  if(LATPARFLAG & n_cov>0) ab <- ab + t(THETA[(dim_irt+3):(dim_irt+2+n_cov)])%*%X

  obj <- cpp_pGrade(
    GRADE=GRADE,
    EXAM=EXAM,
    THETA_IRT=THETA[1:dim_irt],
    THETA_LAT=THETA[(dim_irt+1):(dim_irt+dim_lat)],
    COVARIATES = X,
    N_GRADES=N_GRADES,
    N_EXAMS=N_EXAMS,
    ABILITY=ab,
    LOGFLAG=LOGFLAG,
    LATPARFLAG=LATPARFLAG)

  if(OUT=="prob"){
    return(obj$prob)
  }else{
    return(obj$gr[1:(dim_irt+dim_lat)])
  }
}
#### test grades probabilities output value and gradient ####
set.seed(123)
abilities <- sort(rnorm(3, 0, 2), decreasing = T)

# without structural parameters
for (abi_index in 1:length(abilities)) {
  for (grade in 1:n_grades) {
    for (exam in 1:n_exams) {
      probGG <- RpGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = FALSE,
        LATPARFLAG=FALSE
      )
      lprobG <- RpGrade(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = TRUE,
        LATPARFLAG=FALSE
      )

      lprobGG <- RpGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = TRUE,
        LATPARFLAG=FALSE
      )

      test_that("pGreaterGrades() val", {
        expect_equal(
          probGG,
          exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade])/(1+exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade]))
        )
      })

      test_that("pGreaterGrades() log val", {
        expect_equal(
          lprobGG,
          log(exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade])/(1+exp(mat[exam,n_grades+1]*abilities[abi_index]-mat[exam,grade])))
        )
      })

      if(grade>1){
        probprevGG <- RpGreaterGrades(
          GRADE = grade-1,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=FALSE
        )

        test_that(paste0("pGreaterGrades() decreases with higher grades | exam ", exam, ", grade ", grade, ", abi_index", abi_index ),{
          expect_true(probprevGG>probGG)
        })



      }

      if(grade < n_grades){
        probnextGG <- RpGreaterGrades(
          GRADE = grade+1,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=FALSE
        )
        tmp <- log(probGG-probnextGG)
        test_that("pGrade() log",{
          expect_equal(lprobG, tmp)
        })
      }else{
        test_that("pGrade() log",{
          expect_equal(lprobG, lprobGG)
        })
      }

      if(abi_index>1){
        probprevGG <- RpGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index-1],
          LOGFLAG = FALSE,
          LATPARFLAG=FALSE
        )

        # test_that(paste0("pGreaterGrades() decreases with lower ability | exam ", exam, ", grade ", grade, ", abi_index", abi_index ),{
        #   expect_true(probprevGG>probGG)
        # })
      }






      FUNprobGG <-function(PAR){
        RpGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA = PAR,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=FALSE
        )
      }
      grGG <- RpGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        OUT="gr",
        LOGFLAG = FALSE,
        LATPARFLAG=FALSE
      )
      test_that("pGreaterGrades() gradient",{
        skip_if_not_installed("numDeriv")
        Rval <- numDeriv::grad(func = FUNprobGG, x = theta[1:(dim_irt+dim_lat)])

        expect_equal(Rval, grGG)
      })

      FUNprobG <-function(PAR){
        RpGrade(
          GRADE = grade,
          EXAM = exam-1,
          THETA = PAR,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=FALSE
        )
      }
      grG <- RpGrade(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        OUT="gr",
        LOGFLAG = FALSE,
        LATPARFLAG=FALSE
      )
      test_that("pGrade() gradient",{
        skip_if_not_installed("numDeriv")
        Rval <- numDeriv::grad(func = FUNprobG, x = theta[1:(dim_irt+dim_lat)])

        expect_equal(Rval, grG)
      })
    }
  }
}

# with structural parameters
for (abi_index in 1:length(abilities)) {
  for (grade in 1:n_grades) {
    for (exam in 1:n_exams) {

      probGG <- RpGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = FALSE,
        LATPARFLAG=TRUE
      )
      lprobG <- RpGrade(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = TRUE,
        LATPARFLAG=TRUE
      )

      lprobGG <- RpGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        LOGFLAG = TRUE,
        LATPARFLAG=TRUE
      )



      if(grade>1){
        probprevGG <- RpGreaterGrades(
          GRADE = grade-1,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=TRUE
        )

        test_that(paste0("pGreaterGrades() decreases with higher grades | exam ", exam, ", grade ", grade, ", abi_index", abi_index ),{
          expect_true(probprevGG>probGG)
        })



      }

      if(grade < n_grades){
        probnextGG <- RpGreaterGrades(
          GRADE = grade+1,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=TRUE
        )
        tmp <- log(probGG-probnextGG)
        test_that("pGrade() log",{
          expect_equal(lprobG, tmp)
        })
      }else{
        test_that("pGrade() log",{
          expect_equal(lprobG, lprobGG)
        })
      }

      if(abi_index>1){
        probprevGG <- RpGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=TRUE
        )

        # test_that(paste0("pGreaterGrades() decreases with lower ability | exam ", exam, ", grade ", grade, ", abi_index", abi_index ),{
        #   expect_true(probprevGG>probGG)
        # })
      }






      FUNprobGG <-function(PAR){
        RpGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA = PAR,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=TRUE
        )
      }
      grGG <- RpGreaterGrades(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        OUT="gr",
        LOGFLAG = FALSE,
        LATPARFLAG=TRUE
      )
      test_that("pGreaterGrades() gradient",{
        skip_if_not_installed("numDeriv")
        Rval <- numDeriv::grad(func = FUNprobGG, x = theta[1:(dim_irt+dim_lat)])

        expect_equal(Rval, grGG)
      })

      FUNprobG <-function(PAR){
        RpGrade(
          GRADE = grade,
          EXAM = exam-1,
          THETA = PAR,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=TRUE
        )
      }
      grG <- RpGrade(
        GRADE = grade,
        EXAM = exam-1,
        THETA = theta,
        N_GRADES = n_grades,
        N_EXAMS = n_exams,
        ABILITY = abilities[abi_index],
        OUT="gr",
        LOGFLAG = FALSE,
        LATPARFLAG=TRUE
      )
      test_that("pGrade() gradient",{
        skip_if_not_installed("numDeriv")
        Rval <- numDeriv::grad(func = FUNprobG, x = theta[1:(dim_irt+dim_lat)])

        expect_equal(Rval, grG)
      })
    }
  }
}



test_that("pGreaterGrades() log extreme ability", {

  abilities <- c(1000, -1000)
  for (abi_index in 1:length(abilities)) {
    for (grade in 1:n_grades) {
      for (exam in 1:n_exams) {
        prob <- RpGreaterGrades(
          GRADE = grade,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = TRUE,
          LATPARFLAG=FALSE
        )

        # check for finite values
        expect_true(is.finite(prob))

      }
    }
  }


})


test_that("check pGrade() probability space", {

  set.seed(321)
  abilities <- sort(rnorm(3, 0, 2), decreasing = T)
  for (abi_index in 1:length(abilities)) {
    for (exam in 1:n_exams) {
      probs <- rep(NA, n_grades+1)
      for (grade in 0:n_grades) {


        probs[grade+1] <- RpGrade(
          GRADE = grade,
          EXAM = exam-1,
          THETA = theta,
          N_GRADES = n_grades,
          N_EXAMS = n_exams,
          ABILITY = abilities[abi_index],
          LOGFLAG = FALSE,
          LATPARFLAG=FALSE
        )




      }
      expect_equal(sum(probs),1)
    }

  }

})
