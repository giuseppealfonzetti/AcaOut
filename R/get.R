
#' Compute Marginal Loglik
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#' @param GRID Grid of quadrature points to be used.
#' @param WEIGHTS Weights for quadrature points.
#' @param INIT If TRUE, it evaluates the marginal log likelihood at the initial parameter values
#'
#' @export
get_loglik <- function(FIT, GRID=NULL, WEIGHTS=NULL, INIT=FALSE){
  if(is.null(GRID)) GRID <- FIT$grid
  if(is.null(WEIGHTS)) WEIGHTS <- FIT$weights

  if(!(FIT$mod%in%c("full", "grtc"))){
    stop("Model not available. Provide fit object for `full` or `grtc` models.")
  }

  pars <- FIT$fit$par
  if(INIT){
    pars <- FIT$start_par
  }

  out <- NA
  if(FIT$mod=="full"){

    out <- CRGRTCM_GH(
      THETA=pars,
      EXAMS_GRADES = FIT$data$gradesMat,
      EXAMS_DAYS = FIT$data$timeMat,
      EXAMS_SET = FIT$data$todoMat,
      EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
      MAX_DAY=FIT$data$fulldata$max_time,
      OUTCOME=FIT$data$outcome,
      EXT_COVARIATES=FIT$data$X,
      YEAR_FIRST=FIT$data$first_year,
      YEAR_LAST=FIT$data$last_year,
      YEAR_LAST_EXAM=FIT$data$yle,
      GRID=GRID,
      WEIGHTS=WEIGHTS,
      YB=FIT$data$yb,
      N_GRADES=FIT$data$n_grades,
      N_EXAMS=FIT$data$n_exams,
      GRFLAG = FALSE,
      ROTGRID = TRUE
    )$ll
  }else if(FIT$mod=="grtc"){
    out <- GRTCM_GH(
      THETA = pars,
      EXAMS_GRADES = FIT$data$gradesMat,
      EXAMS_DAYS = FIT$data$timeMat,
      EXAMS_SET = FIT$data$todoMat,
      EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
      MAX_DAY=FIT$data$fulldata$max_time,
      GRID=GRID,
      WEIGHTS=WEIGHTS,
      N_GRADES = FIT$data$n_grades,
      N_EXAMS = FIT$data$n_exams,
      GRFLAG = FALSE,
      ROTGRID = TRUE
    )$ll
  }

  return(out)
}


#' get fitted outcome probabilities given exam performances
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#'
#' @export
get_fitted_outcome <- function(FIT){
  expand_grid(subject_id=FIT$data$obs_ids, year = 1:FIT$data$yb) |>
    mutate(probs=map2(subject_id, year,~{
      ID <- which(FIT$data$obs_ids==.x)
      if(.y < FIT$data$first_year[ID]){
        pout <- rep(NA,4)
        names(pout) <- c("enrollment", "graduation", "dropout", "transfer")
        return(as_tibble(t(pout)))
      }else{
        pout <- sapply(0:3, function(outcome){
          exp(CRGRTCM_GH2(
            THETA = FIT$fit$par,
            EXAMS_GRADES = t(as.matrix(FIT$data$gradesMat[ID,])),
            EXAMS_DAYS = t(as.matrix(FIT$data$timeMat[ID,])),
            EXAMS_SET = t(as.matrix(FIT$data$todoMat[ID,])),
            EXAMS_OBSFLAG = t(as.matrix(!is.na(FIT$data$timeMat[ID,]))),
            MAX_DAY = as.vector(FIT$data$fulldata$max_time[ID]),
            OUTCOME = as.vector(outcome),
            EXT_COVARIATES = t(as.matrix(FIT$data$X[ID,])),
            YEAR_FIRST = as.vector(FIT$data$first_year[ID]),
            YEAR_LAST = as.vector(.y),
            YEAR_LAST_EXAM = as.vector(FIT$data$yle[ID]),
            GRID = FIT$grid,
            WEIGHTS = FIT$weights,
            YB = FIT$data$yb,
            N_GRADES = FIT$data$n_grades,
            N_EXAMS = FIT$data$n_exams,
            GRFLAG = FALSE,
            ROTGRID = TRUE
          )$ll)
        })
        pout <- pout/sum(pout)
        names(pout) <- c("enrollment", "graduation", "dropout", "transfer")
        return(as_tibble(t(pout)))
      }}))|> unnest("probs")
  }
