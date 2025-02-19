#' Compute standard errors from GRTCM fit
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#' @param TIDY Return tidy parameter table
#'
#' @importFrom rlang .data
#' @importFrom dplyr mutate
#'
#' @export
compute_stderr <- function(FIT, NUMDER=FALSE, TIDY = TRUE, GRID=NULL, WEIGHTS=NULL){

  if(!(FIT$mod%in%c("full", "grtc", "ccr"))){
    stop("Model not available. Provide fit object for `full` , `grtc` or `ccr` models.")
  }

  out <- list()
  internal_invH <- FIT$fit$invhessian
  internal_grid <- GRID
  internal_weights <- WEIGHTS
  if(is.null(GRID)) internal_grid <- FIT$grid
  if(is.null(WEIGHTS)) internal_weights <- FIT$weights

  message(paste0("Computing standard errors of ", FIT$mod, " model."))

  if(is.null(FIT$fit$invhessian)) {
    message("Inverse hessian approximation not found in the FIT object. Proceeding with numerical derivatives...")
    NUMDER=TRUE
  }

  dim_irt <- FIT$data$n_exams * (FIT$data$n_grades+3)
  dim_cr <- (FIT$data$yb + FIT$data$cr_ext_cov + 2)*2+1

  if(FIT$mod=="ccr"){
    FIT$fit$par <- c(rep(NA, dim_irt+2), FIT$fit$par)
  }else if(FIT$mod=="grtc"){
    FIT$fit$par <- c(FIT$fit$par, rep(NA, dim_cr))
  }

  reparJacob <- numDeriv::jacobian(func = parVec2Repar, x = FIT$fit$par,
                                   YB=FIT$data$yb,
                                   N_COV=FIT$data$cr_ext_cov,
                                   N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
                                   LABS_EXAMS = FIT$data$exams_labs,
                                   LABS_GRADES = FIT$data$grades_labs,
                                   LABS_COV = colnames(FIT$data$X)
  )

  if(NUMDER){
      NGR <- function(x){
        -cpp_GQ(
          THETA = x,
          EXAMS_GRADES = FIT$data$gradesMat,
          EXAMS_DAYS = FIT$data$timeMat,
          EXAMS_SET = FIT$data$todoMat,
          EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
          COVARIATES = as.matrix(FIT$data$X),
          MAX_DAY = FIT$data$fulldata$max_time,
          OUTCOME = FIT$data$outcome,
          YEAR_FIRST = FIT$data$first_year,
          YEAR_LAST = FIT$data$last_year,
          YEAR_LAST_EXAM = FIT$data$yle,
          YB = FIT$data$yb,
          GRID = internal_grid,
          WEIGHTS = internal_weights,
          N_GRADES = FIT$data$n_grades,
          N_EXAMS = FIT$data$n_exams,
          GRFLAG = TRUE,
          LATPARFLAG = TRUE,
          MOD=FIT$mod
        )$gr}


      numHess <- numDeriv::jacobian(func = NGR, x=FIT$fit$par)
      out[["numHess"]] <- numHess
      internal_invH <- MASS::ginv(numHess)
  }



  out[["invHess"]] <- internal_invH

    if(FIT$mod=="ccr"){
      reparJacob <- reparJacob[-c(1:(dim_irt+2)), -c(1:(dim_irt+2))]
      seVec <- c(rep(NA, dim_irt+2), sqrt(diag(t(reparJacob) %*% internal_invH %*% reparJacob)))
    }else if(FIT$mod=="grtc"){
      reparJacob <- reparJacob[c(1:(dim_irt+2)), c(1:(dim_irt+2))]
      seVec <- c(sqrt(diag(t(reparJacob) %*% internal_invH %*% reparJacob)), rep(NA, dim_cr))
    }else if(FIT$mod=="full"){
      seVec <- sqrt(diag(t(reparJacob) %*% internal_invH %*% reparJacob))
    }



  if(TIDY){
    out[["se"]] <- parVec2Repar(FIT$fit$par,YB=FIT$data$yb,
                        N_COV=FIT$data$cr_ext_cov,
                        N_GRADES = FIT$data$n_grades, N_EXAMS = FIT$data$n_exams,
                        LABS_EXAMS = FIT$data$exams_labs,
                        LABS_GRADES = FIT$data$grades_labs,
                        LABS_COV = colnames(FIT$data$X), TIDY = TRUE) |>
      mutate(
        se  = seVec,
        lb = .data$par-1.96*.data$se,
        ub = .data$par+1.96*.data$se,
        sig = !(.data$lb<0&.data$ub>0)
      )
  }else{
    out[["se"]] <- seVec
  }

  return(out)

}
