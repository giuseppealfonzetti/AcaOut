#' Fit CRGRTCM
#'
#' Fit Competing Risk Graded Response Time Censored Model
#'
#' @param DATA Dataset in long form
#' @param GRADE_COL Label of the column where to find the grades
#' @param EXAM_COL Label of the column where to find exams ids
#' @param STUDID_COL Label of the column where to find students ids.
#' @param TIME_COL Label of the column where to find exams times.
#' @param MAXTIME_COL Label of the column where to find the maximum time observed for each student.
#' @param GRID Grid of quadrature points to be used.
#' @param WEIGHTS Weights for quadrature points
#' @param THETA_START Optional starting parameter vector
#' @param ... additional arguments to be passed to [ucminf::ucminf()]
#'
#' @importFrom stats sd
#' @export
fit_BFGS <- function(DATA,
                        GRADE_COL = "xcat",
                        EXAM_COL = "itemID",
                        STUDID_COL = "subjectID",
                        TIME_COL = "time",
                        MAXTIME_COL = "maxTime",
                        GRID, WEIGHTS, THETA_START = NULL,
                        MOD,
                        ...){

  if(!(MOD%in%c("full", "grtc"))){
    stop("Select MOD among `full` or `grtc`.")
  }

  if(is.null(THETA_START)){
    startIRTMat <- matrix(NA, DATA$n_exams, DATA$n_grades+3)
    startIRTMat[,1:DATA$n_grades] <- matrix(rep(seq(-2,4, length.out = DATA$n_grades), DATA$n_exams), DATA$n_exams, DATA$n_grades, byrow = TRUE)
    startIRTMat[,DATA$n_grades+1] <- 1
    startIRTMat[,DATA$n_grades+2] <- apply(DATA$timeMat, MARGIN = 2, FUN = function(x) mean(log(x), na.rm = TRUE))
    startIRTMat[,DATA$n_grades+3] <- apply(DATA$timeMat, MARGIN = 2, FUN = function(x) 1/sd(x, na.rm = TRUE))
    startLatMat <- diag(1, 2, 2)
    startBeta <- matrix(-1, DATA$yb+DATA$cr_ext_cov+2, 2)
    startGradInt <- 3
    start_par <- parList2Vec(list("irt"=startIRTMat, 'lat_var'=startLatMat, "cr"=list("beta"=startBeta, "grad"=startGradInt)))
  }else{
    start_par <- THETA_START
  }

  NLL <- function(x){
    -cpp_GQ(
      THETA = x,
      EXAMS_GRADES = DATA$gradesMat,
      EXAMS_DAYS = DATA$timeMat,
      EXAMS_SET = DATA$todoMat,
      EXAMS_OBSFLAG = !is.na(DATA$timeMat),
      COVARIATES = as.matrix(DATA$X),
      MAX_DAY = DATA$fulldata$max_time,
      OUTCOME = DATA$outcome,
      YEAR_FIRST = DATA$first_year,
      YEAR_LAST = DATA$last_year,
      YEAR_LAST_EXAM = DATA$yle,
      YB = DATA$yb,
      GRID = GRID,
      WEIGHTS = WEIGHTS,
      N_GRADES = DATA$n_grades,
      N_EXAMS = DATA$n_exams,
      GRFLAG = FALSE,
      LATPARFLAG = TRUE,
      MOD=MOD
    )$ll
  }
  # NLL <- function(x){
  #   -CRGRTCM_GH(
  #     THETA = x,
  #     EXAMS_GRADES = DATA$gradesMat,
  #     EXAMS_DAYS = DATA$timeMat,
  #     EXAMS_SET = DATA$todoMat,
  #     EXAMS_OBSFLAG = !is.na(DATA$timeMat),
  #     MAX_DAY = DATA$fulldata$max_time,
  #     OUTCOME = DATA$outcome,
  #     EXT_COVARIATES = as.matrix(DATA$X) ,
  #     YEAR_FIRST = DATA$first_year,
  #     YEAR_LAST = DATA$last_year,
  #     YEAR_LAST_EXAM = DATA$yle,
  #     GRID = GRID,
  #     WEIGHTS = WEIGHTS,
  #     YB = DATA$yb,
  #     N_GRADES = DATA$n_grades,
  #     N_EXAMS = DATA$n_exams,
  #     GRFLAG = FALSE,
  #     ROTGRID = TRUE
  #   )$ll
  # }

  # NGR <- function(x){
  #   -CRGRTCM_GH(
  #     THETA = x,
  #     EXAMS_GRADES = DATA$gradesMat,
  #     EXAMS_DAYS = DATA$timeMat,
  #     EXAMS_SET = DATA$todoMat,
  #     EXAMS_OBSFLAG = !is.na(DATA$timeMat),
  #     MAX_DAY = DATA$fulldata$max_time,
  #     OUTCOME = DATA$outcome,
  #     EXT_COVARIATES = as.matrix(DATA$X) ,
  #     YEAR_FIRST = DATA$first_year,
  #     YEAR_LAST = DATA$last_year,
  #     YEAR_LAST_EXAM = DATA$yle,
  #     GRID = GRID,
  #     WEIGHTS = WEIGHTS,
  #     YB = DATA$yb,
  #     N_GRADES = DATA$n_grades,
  #     N_EXAMS = DATA$n_exams,
  #     GRFLAG = TRUE,
  #     ROTGRID = TRUE
  #   )$gr
  # }
  NGR <- function(x){
    -cpp_GQ(
      THETA = x,
      EXAMS_GRADES = DATA$gradesMat,
      EXAMS_DAYS = DATA$timeMat,
      EXAMS_SET = DATA$todoMat,
      EXAMS_OBSFLAG = !is.na(DATA$timeMat),
      COVARIATES = as.matrix(DATA$X),
      MAX_DAY = DATA$fulldata$max_time,
      OUTCOME = DATA$outcome,
      YEAR_FIRST = DATA$first_year,
      YEAR_LAST = DATA$last_year,
      YEAR_LAST_EXAM = DATA$yle,
      YB = DATA$yb,
      GRID = GRID,
      WEIGHTS = WEIGHTS,
      N_GRADES = DATA$n_grades,
      N_EXAMS = DATA$n_exams,
      GRFLAG = TRUE,
      LATPARFLAG = TRUE,
      MOD=MOD
    )$gr
  }
  fit <- ucminf::ucminf(par = start_par, fn = NLL, gr = NGR, hessian = 2, ...)

  message(fit$message)

  return(
    list(
      "grid" = GRID,
      "weights" = WEIGHTS,
      "data" = DATA,
      "start_par" = start_par,
      "mod" = MOD,
      "fit" = fit
    ))
}
