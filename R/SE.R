#' Compute standard errors from GRTCM fit
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#' @param METHOD Choose among `sample` `num` and `bfgs`.
#' @param TIDY Return tidy parameter table
#' @param GRID GH points grid
#' @param WEIGHTS GH weights
#' @param VERBOSE TRUE for verbose output
#'
#' @importFrom rlang .data
#' @importFrom dplyr mutate
#'
#' @export
compute_stderr <- function(
  FIT,
  METHOD = c("sample", "num", "bfgs"),
  TIDY = TRUE,
  GRID = NULL,
  WEIGHTS = NULL,
  VERBOSE = FALSE
) {
  if (!(FIT$mod %in% c("full", "grtc", "ccr"))) {
    stop(
      "Model not available. Provide fit object for `full` , `grtc` or `ccr` models."
    )
  }

  out <- list()
  internal_invH <- FIT$fit$invhessian
  internal_grid <- GRID
  internal_weights <- WEIGHTS
  if (is.null(GRID)) {
    internal_grid <- FIT$grid
  }
  if (is.null(WEIGHTS)) {
    internal_weights <- FIT$weights
  }

  if (VERBOSE) {
    message(paste0("Computing standard errors of ", FIT$mod, " model."))
  }

  dim_irt <- FIT$data$par_dims$grtcm
  dim_lat <- FIT$data$par_dims$lat
  dim_cr <- FIT$data$par_dims$cr

  if (FIT$mod == "ccr") {
    FIT$fit$par <- c(rep(NA, dim_irt + dim_lat), FIT$fit$par)
  } else if (FIT$mod == "grtc") {
    FIT$fit$par <- c(FIT$fit$par[1:(dim_irt + dim_lat)], rep(NA, dim_cr))
  }

  reparJacob <- numDeriv::jacobian(
    func = parVec2Repar,
    x = FIT$fit$par,
    YB = FIT$data$data_dims$yb,
    N_COV = FIT$data$data_dims$n_cov,
    N_GRADES = FIT$data$data_dims$n_grades,
    N_EXAMS = FIT$data$data_dims$n_exams,
    LABS_EXAMS = FIT$data$labs$exams,
    LABS_GRADES = FIT$data$labs$grades,
    LABS_COV = FIT$data$labs$cov
  )

  out[["reparJacob"]] <- reparJacob

  if (METHOD == "bfgs") {
    if (is.null(FIT$fit$invhessian)) {
      message(
        "Inverse hessian approximation not found in the FIT object. Proceeding with sample estiator..."
      )
      METHOD <- "sample"
    }
  }

  if (METHOD == "sample") {
    H <- as.matrix(
      Matrix::nearPD(
        cpp_GQ(
          THETA = FIT$fit$par,
          EXAMS_GRADES = FIT$data$gradesMat,
          EXAMS_DAYS = FIT$data$timeMat,
          EXAMS_SET = FIT$data$todoMat,
          EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
          COVARIATES = as.matrix(FIT$data$X),
          MAX_DAY = FIT$data$max_time,
          OUTCOME = FIT$data$outcome,
          YEAR_FIRST = FIT$data$first_year,
          YEAR_LAST = FIT$data$last_year,
          YEAR_LAST_EXAM = FIT$data$yle,
          YB = FIT$data$data_dims$yb,
          GRID = internal_grid,
          WEIGHTS = internal_weights,
          N_GRADES = FIT$data$data_dims$n_grades,
          N_EXAMS = FIT$data$data_dims$n_exams,
          GRFLAG = TRUE,
          LATPARFLAG = TRUE,
          MOD = FIT$mod,
          HFLAG = TRUE
        )$H
      )$mat
    )
    out[["sampleHess"]] <- H
    internal_invH <- as.matrix(Matrix::nearPD(MASS::ginv(H))$mat)
  } else if (METHOD == "num") {
    NGR <- function(x) {
      -cpp_GQ(
        THETA = x,
        EXAMS_GRADES = FIT$data$gradesMat,
        EXAMS_DAYS = FIT$data$timeMat,
        EXAMS_SET = FIT$data$todoMat,
        EXAMS_OBSFLAG = !is.na(FIT$data$timeMat),
        COVARIATES = as.matrix(FIT$data$X),
        MAX_DAY = FIT$data$max_time,
        OUTCOME = FIT$data$outcome,
        YEAR_FIRST = FIT$data$first_year,
        YEAR_LAST = FIT$data$last_year,
        YEAR_LAST_EXAM = FIT$data$yle,
        YB = FIT$data$data_dims$yb,
        GRID = internal_grid,
        WEIGHTS = internal_weights,
        N_GRADES = FIT$data$data_dims$n_grades,
        N_EXAMS = FIT$data$data_dims$n_exams,
        GRFLAG = TRUE,
        LATPARFLAG = TRUE,
        MOD = FIT$mod
      )$gr
    }

    numHess <- numDeriv::jacobian(func = NGR, x = FIT$fit$par)
    out[["numHess"]] <- numHess
    internal_invH <- MASS::ginv(numHess)
  }

  out[["invHess"]] <- internal_invH

  if (FIT$mod == "ccr") {
    reparJacob <- reparJacob[-c(1:(dim_irt + 2)), -c(1:(dim_irt + 2))]
    seVec <- c(
      rep(NA, dim_irt + 2),
      sqrt(diag(t(reparJacob) %*% internal_invH %*% reparJacob))
    )
  } else if (FIT$mod == "grtc") {
    reparJacob <- reparJacob[c(1:(dim_irt + 2)), c(1:(dim_irt + 2))]
    seVec <- c(
      sqrt(diag(t(reparJacob) %*% internal_invH %*% reparJacob)),
      rep(NA, dim_cr)
    )
  } else if (FIT$mod == "full") {
    seVec <- sqrt(diag(t(reparJacob) %*% internal_invH %*% reparJacob))
    seVec[!is.finite(seVec)] <- 1
  }

  if (TIDY) {
    out[["se"]] <- parVec2Repar(
      FIT$fit$par,
      YB = FIT$data$data_dims$yb,
      N_COV = FIT$data$data_dims$n_cov,
      N_GRADES = FIT$data$data_dims$n_grades,
      N_EXAMS = FIT$data$data_dims$n_exams,
      LABS_EXAMS = FIT$data$labs$exams,
      LABS_GRADES = FIT$data$labs$grades,
      LABS_COV = FIT$data$labs$cov,
      TIDY = TRUE
    ) |>
      mutate(
        se = seVec,
        lb = .data$par - 1.96 * .data$se,
        ub = .data$par + 1.96 * .data$se,
        sig = !(.data$lb < 0 & .data$ub > 0)
      )
  } else {
    out[["se"]] <- seVec
  }

  return(out)
}
