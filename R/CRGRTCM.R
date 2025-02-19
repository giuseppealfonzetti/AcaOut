#' #' Fit CRGRTCM
#' #'
#' #' Fit Competing Risk Graded Response Time Censored Model
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
#' #' @param ... additional arguments to be passed to [ucminf::ucminf()]
#' #'
#' #' @importFrom stats sd
#' #' @export
#' fit_CRGRTCM <- function(DATA,
#'                         GRADE_COL = "xcat",
#'                         EXAM_COL = "itemID",
#'                         STUDID_COL = "subjectID",
#'                         TIME_COL = "time",
#'                         MAXTIME_COL = "maxTime",
#'                         GRID, WEIGHTS, THETA_START = NULL,
#'                         ...){
#'
#'
#'
#'   if(is.null(THETA_START)){
#'     startIRTMat <- matrix(NA, DATA$n_exams, DATA$n_grades+3)
#'     startIRTMat[,1:DATA$n_grades] <- matrix(rep(seq(-2,4, length.out = DATA$n_grades), DATA$n_exams), DATA$n_exams, DATA$n_grades, byrow = TRUE)
#'     startIRTMat[,DATA$n_grades+1] <- 1
#'     startIRTMat[,DATA$n_grades+2] <- apply(DATA$timeMat, MARGIN = 2, FUN = function(x) mean(log(x), na.rm = TRUE))
#'     startIRTMat[,DATA$n_grades+3] <- apply(DATA$timeMat, MARGIN = 2, FUN = function(x) 1/sd(x, na.rm = TRUE))
#'     startLatMat <- diag(1, 2, 2)
#'     startBeta <- matrix(rnorm(2*(DATA$yb+DATA$cr_ext_cov+2), 0, .1), DATA$yb+DATA$cr_ext_cov+2, 2)
#'     startGradInt <- 3
#'     start_par <- parList2Vec(list("irt"=startIRTMat, 'lat_var'=startLatMat, "cr"=list("beta"=startBeta, "grad"=startGradInt)))
#'   }else{
#'     start_par <- THETA_START
#'   }
#'
#'   NLL <- function(x){
#'     -CRGRTCM_GH(
#'       THETA = x,
#'       EXAMS_GRADES = DATA$gradesMat,
#'       EXAMS_DAYS = DATA$timeMat,
#'       EXAMS_SET = DATA$todoMat,
#'       EXAMS_OBSFLAG = !is.na(DATA$timeMat),
#'       MAX_DAY = DATA$fulldata$max_time,
#'       OUTCOME = DATA$outcome,
#'       EXT_COVARIATES = as.matrix(DATA$X) ,
#'       YEAR_FIRST = DATA$first_year,
#'       YEAR_LAST = DATA$last_year,
#'       YEAR_LAST_EXAM = DATA$yle,
#'       GRID = GRID,
#'       WEIGHTS = WEIGHTS,
#'       YB = DATA$yb,
#'       N_GRADES = DATA$n_grades,
#'       N_EXAMS = DATA$n_exams,
#'       GRFLAG = FALSE,
#'       ROTGRID = TRUE
#'     )$ll
#'   }
#'
#'   NGR <- function(x){
#'     -CRGRTCM_GH(
#'       THETA = x,
#'       EXAMS_GRADES = DATA$gradesMat,
#'       EXAMS_DAYS = DATA$timeMat,
#'       EXAMS_SET = DATA$todoMat,
#'       EXAMS_OBSFLAG = !is.na(DATA$timeMat),
#'       MAX_DAY = DATA$fulldata$max_time,
#'       OUTCOME = DATA$outcome,
#'       EXT_COVARIATES = as.matrix(DATA$X) ,
#'       YEAR_FIRST = DATA$first_year,
#'       YEAR_LAST = DATA$last_year,
#'       YEAR_LAST_EXAM = DATA$yle,
#'       GRID = GRID,
#'       WEIGHTS = WEIGHTS,
#'       YB = DATA$yb,
#'       N_GRADES = DATA$n_grades,
#'       N_EXAMS = DATA$n_exams,
#'       GRFLAG = TRUE,
#'       ROTGRID = TRUE
#'     )$gr
#'   }
#'
#'   fit <- ucminf::ucminf(par = start_par, fn = NLL, gr = NGR, hessian = 2, ...)
#'
#'   # fit <- optimx::optimx(par=start_par, fn = NLL, gr = NGR,  ...)
#'   return(
#'     list(
#'       "data" = DATA,
#'       "start_par" = start_par,
#'       "mod" = "crgrtcm",
#'       "fit" = fit
#'   ))
#' }
#'
#'
#' #' MAP GRTCM
#' #'
#' #' Compute MAP scores for Graded Response Time Censored Model
#' #'
#' #' @param FIT Output from [fit_CRGRTCM]
#' #' @param TIDY Return tidy parameter table
#' #' @param MATSTART Starting values for the latent scores
#' #'
#' #' @importFrom dplyr as_tibble
#' #' @export
#' map_CRGRTCM <- function(FIT, MATSTART=NULL, TIDY = TRUE){
#'
#'   if(is.null(MATSTART)){
#'     MATSTART <- matrix(0, length(FIT$data$obs_ids), 2)
#'   }
#'
#'
#'   ID_COMPLETENLL <- function(LAT, ID){
#'     -CRGRTCM_complete_obs(
#'       THETA = FIT$fit$par,
#'       EXAMS_GRADES = FIT$data$gradesMat[ID,],
#'       EXAMS_DAYS = FIT$data$timeMat[ID,],
#'       EXAMS_SET = FIT$data$todoMat[ID,],
#'       EXAMS_OBSFLAG = !is.na(FIT$data$timeMat[ID,]),
#'       OUTCOME=FIT$data$outcome[ID],
#'       YEAR_FIRST=FIT$data$first_year[ID],
#'       YEAR_LAST=FIT$data$last_year[ID],
#'       YEAR_LAST_EXAM=FIT$data$yle[ID],
#'       EXT_COVARIATES=FIT$data$X[ID,],
#'       YB=FIT$data$yb,
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
#'                          fit <- ucminf::ucminf(par = as.numeric(MATSTART[ID,]), fn = ID_COMPLETENLL, ID = ID)
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
#' #' Compute standard errors from GRTCM fit
#' #'
#' #' @param FIT Output from [fit_CRGRTCM]
#' #' @param TIDY Return tidy parameter table
#' #'
#' #' @importFrom rlang .data
#' #' @importFrom dplyr mutate
#' #'
#' #' @export
#' compute_se <- function(FIT, TIDY = TRUE){
#'
#'   dim_irt <- FIT$data$n_exams * (FIT$data$n_grades+3)
#'   dim_cr <- (FIT$data$yb + FIT$data$cr_ext_cov + 2)*2+1
#'
#'   if(FIT$mod=="ccr"){
#'     FIT$fit$par <- c(rep(NA, dim_irt+2), FIT$fit$par)
#'   }else if(FIT$mod=="grtc"){
#'     FIT$fit$par <- c(FIT$fit$par, rep(NA, dim_cr))
#'   }
#'
#'   reparJacob <- numDeriv::jacobian(func = parVec2Repar, x = FIT$fit$par,
#'                                    YB=FIT$data$yb,
#'                                    N_COV=FIT$data$cr_ext_cov,
#'                                    N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
#'                                    LABS_EXAMS = FIT$data$exams_labs,
#'                                    LABS_GRADES = FIT$data$grades_labs,
#'                                    LABS_COV = colnames(FIT$data$X)
#'   )
#'
#'   if(FIT$mod=="ccr"){
#'     reparJacob <- reparJacob[-c(1:(dim_irt+2)), -c(1:(dim_irt+2))]
#'     seVec <- c(rep(NA, dim_irt+2), sqrt(diag(t(reparJacob) %*% FIT$fit$invhessian %*% reparJacob)))
#'   }else if(FIT$mod=="grtcm"){
#'     reparJacob <- reparJacob[c(1:(dim_irt+2)), c(1:(dim_irt+2))]
#'     seVec <- c(sqrt(diag(t(reparJacob) %*% FIT$fit$invhessian %*% reparJacob)), rep(NA, dim_cr))
#'   }else{
#'     seVec <- sqrt(diag(t(reparJacob) %*% FIT$fit$invhessian %*% reparJacob))
#'   }
#'
#'   if(TIDY){
#'     out <- parVec2Repar(FIT$fit$par,YB=FIT$data$yb,
#'                         N_COV=FIT$data$cr_ext_cov,
#'                         N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
#'                         LABS_EXAMS = FIT$data$exams_labs,
#'                         LABS_GRADES = FIT$data$grades_labs,
#'                         LABS_COV = colnames(FIT$data$X), TIDY = TRUE) |>
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
