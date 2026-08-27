#' Compute standard errors from GRTCM fit
#'
#' @param FIT Output from [fit_EM] or [fit_BFGS]
#' @param METHOD One of `sample` (OPG), `num` (inverse observed information),
#'   `bfgs` (optimiser inverse-Hessian) or `sandwich` (robust Huber-White
#'   estimator, consistent under misspecification).
#' @param TIDY Return tidy parameter table
#' @param GRID GH points grid
#' @param WEIGHTS GH weights
#' @param NCORES Cores for the finite-difference observed information used by
#'   `METHOD="sandwich"`. Defaults to `detectCores() - 1`, or `1` on Windows.
#' @param VERBOSE TRUE for verbose output
#'
#' @importFrom rlang .data
#' @importFrom dplyr mutate
#'
#' @export
compute_stderr <- function(
  FIT,
  METHOD = c("sample", "num", "bfgs", "sandwich"),
  TIDY = TRUE,
  GRID = NULL,
  WEIGHTS = NULL,
  NCORES = NULL,
  VERBOSE = FALSE
) {
  METHOD <- match.arg(METHOD)
  ncores <- if (is.null(NCORES)) {
    if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)
  } else {
    as.integer(NCORES)
  }
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

  gq_eval <- function(theta, HFLAG = FALSE) {
    cpp_GQ(
      THETA = theta,
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
      HFLAG = HFLAG
    )
  }

  fd_information <- function(theta, h_rel = 1e-4, cores = 1L) {
    p <- length(theta)
    col_j <- function(j) {
      hj <- h_rel * max(abs(theta[j]), 1)
      tp <- theta
      tm <- theta
      tp[j] <- tp[j] + hj
      tm[j] <- tm[j] - hj
      -(gq_eval(tp)$gr - gq_eval(tm)$gr) / (2 * hj)
    }
    cols <- if (cores > 1L) {
      parallel::mclapply(seq_len(p), col_j, mc.cores = cores)
    } else {
      lapply(seq_len(p), col_j)
    }
    A <- do.call(cbind, cols)
    (A + t(A)) / 2
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
    H <- as.matrix(Matrix::nearPD(gq_eval(FIT$fit$par, HFLAG = TRUE)$H)$mat)
    out[["sampleHess"]] <- H
    internal_invH <- as.matrix(Matrix::nearPD(MASS::ginv(H))$mat)
  } else if (METHOD == "num") {
    numHess <- numDeriv::jacobian(
      func = function(x) -gq_eval(x)$gr,
      x = FIT$fit$par
    )
    out[["numHess"]] <- numHess
    internal_invH <- MASS::ginv(numHess)
  } else if (METHOD == "sandwich") {
    B <- gq_eval(FIT$fit$par, HFLAG = TRUE)$H
    B <- (B + t(B)) / 2
    if (any(!is.finite(B))) {
      warning("non-finite entries in the OPG meat B; check for students with underflowing likelihood")
    }
    out[["meat"]] <- B

    A <- fd_information(FIT$fit$par, cores = ncores)
    out[["bread"]] <- A

    Ainv <- MASS::ginv(A)
    internal_invH <- Ainv %*% B %*% Ainv
    internal_invH <- (internal_invH + t(internal_invH)) / 2
    out[["sandwichV"]] <- internal_invH
  }

  out[["invHess"]] <- internal_invH

  if (FIT$mod == "ccr") {
    reparJacob <- reparJacob[-c(1:(dim_irt + 2)), -c(1:(dim_irt + 2))]
  } else if (FIT$mod == "grtc") {
    reparJacob <- reparJacob[c(1:(dim_irt + 2)), c(1:(dim_irt + 2))]
  }
  seRepar <- sqrt(rowSums((reparJacob %*% internal_invH) * reparJacob))

  if (any(!is.finite(seRepar))) {
    warning(
      sum(!is.finite(seRepar)),
      " of ",
      length(seRepar),
      " standard errors are not finite",
      call. = FALSE
    )
  }

  seVec <- switch(
    FIT$mod,
    ccr = c(rep(NA, dim_irt + 2), seRepar),
    grtc = c(seRepar, rep(NA, dim_cr)),
    full = seRepar
  )

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
