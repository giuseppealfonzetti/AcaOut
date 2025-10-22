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
compute_eap <- function(FIT, TIDY = TRUE){

  if(!(FIT$mod%in%c("full", "grtc"))){
    stop("Model not available. Provide fit object for `full` or `grtc` models.")
  }

  message(paste0("Computing EAP latent score estimates of ", FIT$mod, " model."))


  mat <- cpp_EAP(
    THETA = FIT$fit$par,
    EXAMS_GRADES = FIT$data$gradesMat,
    EXAMS_DAYS = FIT$data$timeMat,
    EXAMS_SET = FIT$data$todoMat,
    EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
    MAX_DAY = FIT$data$max_time,
    OUTCOME = FIT$data$outcome,
    EXT_COVARIATES = as.matrix(FIT$data$X) ,
    YEAR_FIRST = FIT$data$first_year,
    YEAR_LAST = FIT$data$last_year,
    YEAR_LAST_EXAM = FIT$data$yle,
    GRID = FIT$grid,
    WEIGHTS = FIT$weights,
    YB = FIT$data$data_dims$yb,
    N_GRADES = FIT$data$data_dims$n_grades,
    N_EXAMS = FIT$data$data_dims$n_exams,
    MOD = FIT$mod,
    VERBOSE=FALSE
    )$EAP





  rownames(mat) <- FIT$data$labs$obs
  colnames(mat) <- c("ability", "speed")
  if(TIDY){
    mat <- as_tibble(mat, rownames = "subject_id")
  }

  return(mat)
}
