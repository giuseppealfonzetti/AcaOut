#' #' Fit GRTCM
#' #'
#' #' Fit Graded Response Time Censored Model
#' #'
#' #' @param DATA Dataset in long form
#' #' @param GRADE_COL Label of the column where to find the grades
#' #' @param EXAM_COL Label of the column where to find exams ids
#' #' @param STUDID_COL Label of the column where to find students ids.
#' #' @param TIME_COL Label of the column where to find exams times.
#' #' @param MAXTIME_COL Label of the column where to find the maximum time observed for each student.
#' #' @param GRID Grid of quadrature points to be used.
#' #' @param WEIGHTS Weights for quadrature points
#' #' @param THETA_START Optional starting parameter vector
#' #'
#' #' @importFrom stats sd
#' #' @export
#' fit_GRTCM <- function(DATA, GRADE_COL = "xcat",
#'                       EXAM_COL = "itemID",
#'                       STUDID_COL = "subjectID",
#'                       TIME_COL = "time",
#'                       MAXTIME_COL = "maxTime",
#'                       GRID, WEIGHTS, THETA_START = NULL){
#'
#'   # dataStruct <- dataRestruct(DATA_IRT   = IRTDATA,
#'   #                            DATA_CR    = CRDATA,
#'   #                            GRADE_COL  = GRADE_COL,
#'   #                            EXAM_COL   = EXAM_COL,
#'   #                            STUDID_COL = STUDID_COL,
#'   #                            TIME_COL   = TIME_COL,
#'   #                            MAXTIME_COL= MAXTIME_COL)
#'
#'   dim_irt <- DATA$n_exams * (DATA$n_grades+3)
#'   dim_cr <- (DATA$yb + DATA$cr_ext_cov + 2)*2+1
#'
#'   if(is.null(THETA_START)){
#'     startIRTMat <- matrix(NA, DATA$n_exams, DATA$n_grades+3)
#'     startIRTMat[,1:DATA$n_grades] <- matrix(rep(seq(-2,4, length.out = DATA$n_grades), DATA$n_exams), DATA$n_exams, DATA$n_grades, byrow = TRUE)
#'     startIRTMat[,DATA$n_grades+1] <- 1
#'     startIRTMat[,DATA$n_grades+2] <- apply(DATA$timeMat, MARGIN = 2, FUN = function(x) mean(log(x), na.rm = TRUE))
#'     startIRTMat[,DATA$n_grades+3] <- apply(DATA$timeMat, MARGIN = 2, FUN = function(x) 1/sd(x, na.rm = TRUE))
#'     startLatMat <- diag(1, 2, 2)
#'     startBeta <- matrix(0, DATA$yb+DATA$cr_ext_cov+2, 2)
#'     startGradInt <- 0
#'     start_par <- parList2Vec(list("irt"=startIRTMat, 'lat_var'=startLatMat, "cr"=list("beta"=startBeta, "grad"=startGradInt)))
#'   }else{
#'     start_par <- THETA_START
#'   }
#'
#'   NLL <- function(x){
#'     -GRTCM_GH(
#'       THETA = x,
#'       EXAMS_GRADES = DATA$gradesMat,
#'       EXAMS_DAYS = DATA$timeMat,
#'       EXAMS_SET = DATA$todoMat,
#'       EXAMS_OBSFLAG = !is.na(DATA$timeMat),
#'       MAX_DAY = DATA$fulldata$max_time,
#'       GRID = GRID,
#'       WEIGHTS = WEIGHTS,
#'       N_GRADES = DATA$n_grades,
#'       N_EXAMS = DATA$n_exams,
#'       GRFLAG = FALSE,
#'       ROTGRID = TRUE
#'     )$ll
#'   }
#'
#'   NGR <- function(x){
#'     -GRTCM_GH(
#'       THETA = x,
#'       EXAMS_GRADES = DATA$gradesMat,
#'       EXAMS_DAYS = DATA$timeMat,
#'       EXAMS_SET = DATA$todoMat,
#'       EXAMS_OBSFLAG = !is.na(DATA$timeMat),
#'       MAX_DAY = DATA$fulldata$max_time,
#'       GRID = GRID,
#'       WEIGHTS = WEIGHTS,
#'       N_GRADES = DATA$n_grades,
#'       N_EXAMS = DATA$n_exams,
#'       GRFLAG = TRUE,
#'       ROTGRID = TRUE
#'     )$gr[1:(dim_irt+2)]
#'   }
#'
#'   fit <- ucminf::ucminf(par = start_par[1:(dim_irt+2)], fn = NLL, gr = NGR, hessian = 2)
#'
#'   return(
#'     list(
#'       "data" = DATA,
#'       "mod" = "grtcm",
#'       "start_par" = start_par,
#'       "fit" = fit
#'   ))
#' }
#'
#' #' MAP GRTCM
#' #'
#' #' Compute MAP scores for Graded Response Time Censored Model
#' #'
#' #' @param FIT Output from [fit_GRTCM]
#' #' @param TIDY Return tidy parameter table
#' #'
#' #' @importFrom dplyr as_tibble
#' #' @export
#' map_GRTCM <- function(FIT, TIDY = TRUE){
#'
#'
#'
#'   ID_COMPLETENLL <- function(LAT, ID){
#'     -GRTCM_complete_obs(
#'       THETA = FIT$fit$par,
#'       EXAMS_GRADES = FIT$data$gradesMat[ID,],
#'       EXAMS_DAYS = FIT$data$timeMat[ID,],
#'       EXAMS_SET = FIT$data$todoMat[ID,],
#'       EXAMS_OBSFLAG = !is.na(FIT$data$timeMat[ID,]),
#'       MAX_DAY = FIT$data$fulldata$max_time[ID],
#'       N_GRADES = FIT$data$n_grades,
#'       N_EXAMS = FIT$data$n_exams,
#'       ABILITY = LAT[1],
#'       SPEED = LAT[2])
#'   }
#'
#'   mat <- Reduce(rbind,
#'                 lapply(1:nrow(FIT$data$gradesMat),
#'                        function(ID){
#'                          fit <- ucminf::ucminf(par = c(0,0), fn = ID_COMPLETENLL, ID = ID)
#'                          as.numeric(fit$par)}))
#'   rownames(mat) <- FIT$data$obs_ids
#'   colnames(mat) <- c("ability", "speed")
#'   if(TIDY){
#'     mat <- as_tibble(mat, rownames = "subject_id")
#'   }
#'
#'   return(mat)
#' }
#'
#'
#'
#' #' Compute standard errors from GRTCM fit
#' #'
#' #' @param FIT Output from [fit_GRTCM]
#' #' @param TIDY Return tidy parameter table
#' #'
#' #' @importFrom rlang .data
#' #' @importFrom dplyr mutate
#' #'
#' #' @export
#' se_GRTCM <- function(FIT, TIDY = TRUE){
#'
#'   reparJacob <- numDeriv::jacobian(func = parVec2Repar, x = FIT$fit$par,
#'                                    N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
#'                                    LABS_EXAMS = FIT$data$exams_labs, LABS_GRADES = FIT$data$grades_labs)
#'
#'   seVec <- sqrt(diag(t(reparJacob) %*% FIT$fit$invhessian %*% reparJacob))
#'
#'   if(TIDY){
#'     out <- parVec2Repar(FIT$fit$par,
#'                            N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
#'                            LABS_EXAMS = FIT$data$exams_labs, LABS_GRADES = FIT$data$grades_labs, TIDY = TRUE) |>
#'       mutate(
#'         se  = seVec,
#'         lb = .data$par-1.96*.data$se,
#'         ub = .data$par+1.96*.data$se,
#'         sig = !(.data$lb<0&.data$ub>0)
#'       )
#'   }else{
#'     out <- seVec
#'   }
#'
#'   return(out)
#'
#' }
