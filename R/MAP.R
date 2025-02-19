#' MAP GRTCM
#'
#' Compute MAP scores for Graded Response Time Censored Model
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#' @param TIDY Return tidy parameter table
#' @param MATSTART Starting values for the latent scores
#'
#' @importFrom dplyr as_tibble
#' @export
compute_map <- function(FIT, MATSTART=NULL, TIDY = TRUE){

  if(!(FIT$mod%in%c("full", "grtc"))){
    stop("Model not available. Provide fit object for `full` or `grtc` models.")
  }

  message(paste0("Computing MAP latent score estimates of ", FIT$mod, " model."))


  if(is.null(MATSTART)){
    MATSTART <- matrix(0, length(FIT$data$obs_ids), 2)
  }

  if(FIT$mod=="full"){
    ID_COMPLETENLL <- function(LAT, ID){
      -CRGRTCM_complete_obs(
        THETA = FIT$fit$par,
        EXAMS_GRADES = FIT$data$gradesMat[ID,],
        EXAMS_DAYS = FIT$data$timeMat[ID,],
        EXAMS_SET = FIT$data$todoMat[ID,],
        EXAMS_OBSFLAG = !is.na(FIT$data$timeMat[ID,]),
        OUTCOME=FIT$data$outcome[ID],
        YEAR_FIRST=FIT$data$first_year[ID],
        YEAR_LAST=FIT$data$last_year[ID],
        YEAR_LAST_EXAM=FIT$data$yle[ID],
        COVARIATES=FIT$data$X[ID,],
        YB=FIT$data$yb,
        MAX_DAY = FIT$data$fulldata$max_time[ID],
        N_GRADES = FIT$data$n_grades,
        N_EXAMS = FIT$data$n_exams,
        ABILITY = LAT[1],
        SPEED = LAT[2])
    }
  }else if(FIT$mod=="grtc"){
    ID_COMPLETENLL <- function(LAT, ID){
      -GRTCM_complete_obs(
        THETA = FIT$fit$par,
        EXAMS_GRADES = FIT$data$gradesMat[ID,],
        EXAMS_DAYS = FIT$data$timeMat[ID,],
        EXAMS_SET = FIT$data$todoMat[ID,],
        EXAMS_OBSFLAG = !is.na(FIT$data$timeMat[ID,]),
        COVARIATES = FIT$data$X[ID,],
        MAX_DAY = FIT$data$fulldata$max_time[ID],
        N_GRADES = FIT$data$n_grades,
        N_EXAMS = FIT$data$n_exams,
        ABILITY = LAT[1],
        SPEED = LAT[2])
    }
  }




  mat <- Reduce(rbind,
                lapply(1:nrow(FIT$data$gradesMat),
                       function(ID){
                         fit <- ucminf::ucminf(par = as.numeric(MATSTART[ID,]), fn = ID_COMPLETENLL, ID = ID)
                         as.numeric(fit$par)}))
  rownames(mat) <- FIT$data$obs_ids
  colnames(mat) <- c("ability", "speed")
  if(TIDY){
    mat <- as_tibble(mat, rownames = "subject_id")
  }

  return(mat)
}
