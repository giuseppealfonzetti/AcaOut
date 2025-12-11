#' Fit measurement model and joint model with BFGS
#'
#' @param DATA data as per [check_data()].
#' @param GRID Grid of quadrature points to be used.
#' @param WEIGHTS Weights for quadrature points
#' @param THETA_START Optional starting parameter vector
#' @param MOD Choose `full` for joint model or `grtcm` for measurement model.
#' @param M_MAX_ITER Max number of iterations M-STEP.
#' @param MAX_ITER Maximum number of EM iterations.
#' @param TOL tolerance level Q-function improvement.
#' @param VERBOSE TRUE for verbose output.
#' @param THETA_START Optional starting parameter vector.
#'
#' @importFrom stats sd
#' @export
fit_EM <- function(
  DATA,
  GRID,
  WEIGHTS,
  MOD,
  M_MAX_ITER,
  MAX_ITER,
  TOL,
  VERBOSE = TRUE,
  THETA_START = NULL
) {
  if (!(MOD %in% c("full", "grtc"))) {
    stop("Select MOD among `full` or `grtc`.")
  }

  if (is.null(THETA_START)) {
    startIRTMat <- matrix(
      NA,
      DATA$data_dims$n_exams,
      DATA$data_dims$n_grades + 3
    )
    startIRTMat[, 1:DATA$data_dims$n_grades] <- matrix(
      rep(
        seq(-2, 4, length.out = DATA$data_dims$n_grades),
        DATA$data_dims$n_exams
      ),
      DATA$data_dims$n_exams,
      DATA$data_dims$n_grades,
      byrow = TRUE
    )
    startIRTMat[, DATA$data_dims$n_grades + 1] <- 1
    startIRTMat[, DATA$data_dims$n_grades + 2] <- apply(
      DATA$timeMat,
      MARGIN = 2,
      FUN = function(x) mean(log(x), na.rm = TRUE)
    )
    startIRTMat[, DATA$data_dims$n_grades + 3] <- apply(
      DATA$timeMat,
      MARGIN = 2,
      FUN = function(x) 1 / sd(log(x), na.rm = TRUE)
    )
    startLatMat <- diag(1, 2, 2)
    startBeta <- matrix(0, DATA$data_dims$yb + DATA$data_dims$n_cov + 2, 2)
    startGradInt <- 5
    start_par <- parList2Vec(list(
      "irt" = startIRTMat,
      'lat_var' = startLatMat,
      "cr" = list("beta" = startBeta, "grad" = startGradInt)
    ))
  } else {
    start_par <- THETA_START
  }

  fit <- cpp_EM(
    THETA_START = start_par,
    EXAMS_GRADES = DATA$gradesMat,
    EXAMS_DAYS = DATA$timeMat,
    EXAMS_SET = DATA$todoMat,
    EXAMS_OBSFLAG = !is.na(DATA$timeMat),
    MAX_DAY = DATA$max_time,
    OUTCOME = DATA$outcome,
    EXT_COVARIATES = as.matrix(DATA$X),
    YEAR_FIRST = DATA$first_year,
    YEAR_LAST = DATA$last_year,
    YEAR_LAST_EXAM = DATA$yle,
    GRID = GRID,
    WEIGHTS = WEIGHTS,
    YB = DATA$data_dims$yb,
    N_GRADES = DATA$data_dims$n_grades,
    N_EXAMS = DATA$data_dims$n_exams,
    M_MAX_ITER = M_MAX_ITER,
    MAX_ITER = MAX_ITER,
    TOL = TOL,
    MOD = MOD,
    VERBOSE = VERBOSE
  )

  return(
    list(
      "grid" = GRID,
      "weights" = WEIGHTS,
      "data" = DATA,
      "start_par" = start_par,
      "mod" = MOD,
      "fit" = fit
    )
  )
}
