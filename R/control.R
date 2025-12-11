check_grades <- function(MAT) {
  stopifnot(is.matrix(MAT))
  stopifnot(all(is.finite(MAT) | is.na(MAT)))
  stopifnot(all(MAT[is.finite(MAT)] > 0))
  gtab <- apply(MAT, MARGIN = 2, FUN = function(x) dim(table(x)))
  if (!all(gtab == gtab[1])) {
    stop(
      "Some grades have not been observed on some exams. Investigate with `apply(MAT, MARGIN=2, FUN=table)`."
    )
  }

  MAT[is.finite(MAT)] <- as.integer(MAT[is.finite(MAT)])
  return(MAT)
}

check_todo <- function(MAT) {
  stopifnot(is.matrix(MAT))
  stopifnot(all(is.logical(MAT)))
  if (any(colSums(MAT) == 0)) {
    stop(
      "Some exams are not included in any study plan.\n Investigate with `colSums(MAT)`."
    )
  }

  return(MAT)
}

check_times <- function(MAT) {
  stopifnot(is.matrix(MAT))
  stopifnot(all(is.finite(MAT) | is.na(MAT)))
  stopifnot(all(MAT[is.finite(MAT)] > 0))
  ttab <- apply(MAT, MARGIN = 2, function(x) length(x[is.finite(x)]))
  if (any(ttab == 0)) {
    stop(
      "Some exams have not been passed by any student yet. Investigate with `apply(MAT, MARGIN = 2, function(x) length(x[is.finite(x)]))`."
    )
  }

  MAT[is.finite(MAT)] <- as.integer(MAT[is.finite(MAT)])

  return(MAT)
}

check_covariates <- function(MAT) {
  stopifnot(is.matrix(MAT))
  stopifnot(all(is.finite(MAT)))
  xtab <- apply(MAT, MARGIN = 2, function(x) all(x == mean(x)))
  if (any(xtab)) {
    stop("One or more covariates collects one single observed value.")
  }

  return(MAT)
}


check_outcome <- function(VEC) {
  stopifnot(is.vector(VEC))
  stopifnot(all(is.finite(VEC)))
  if (!all(VEC >= 0 & VEC <= 3)) {
    stop(
      "Outcomes not coded correctly. Use 0 for enrollment, 1 for graduation, 2 for dropout, 3 for transfer."
    )
  }
  if (dim(table(VEC)) < 3) {
    stop("One or more outcomes have not been observed on the whole dataset.")
  }
  VEC[is.finite(VEC)] <- as.integer(VEC[is.finite(VEC)])
  return(VEC)
}

check_year <- function(VEC) {
  stopifnot(is.vector(VEC))
  stopifnot(all(is.finite(VEC)))
  stopifnot(all(VEC >= 0))
  VEC[is.finite(VEC)] <- as.integer(VEC[is.finite(VEC)])

  return(VEC)
}

check_yearle <- function(VEC) {
  stopifnot(is.vector(VEC))
  stopifnot(all(is.finite(VEC) | is.na(VEC)))
  stopifnot(all(VEC[is.finite(VEC)] >= 0))
  VEC[is.finite(VEC)] <- as.integer(VEC[is.finite(VEC)])
  VEC[!is.finite(VEC)] <- 100L
  return(VEC)
}
check_maxt <- function(VEC) {
  stopifnot(is.vector(VEC))
  stopifnot(all(is.finite(VEC)))
  stopifnot(all(VEC[is.finite(VEC)] > 0))
  VEC[is.finite(VEC)] <- as.integer(VEC[is.finite(VEC)])
  return(VEC)
}


#' Helper function to check data structure
#'
#' @param GRADES Grades matrix
#' @param TIMES Time matrix
#' @param TODO To-do matrix
#' @param OUTCOME Outcome vector
#' @param X COvariate matrix
#' @param FIRST_YEAR Vector containing start years
#' @param LAST_YEAR Vector containing last years
#' @param LAST_EXAM_YEAR Vector containing last exam year
#' @param MAX_TIME Vector containing max observable day
#' @param LABS_EXAMS Exam labels
#' @param LABS_OBS Observation labels
#' @param LABS_GRADES Grades labels
#' @param LABS_COV Covariates labels
#' @param SUBSET Subset of the sample to be returned
#' @param VERBOSE TRUE for verbose output
#'
#' @export
check_data <- function(
  GRADES,
  TIMES,
  TODO,
  OUTCOME,
  X,
  FIRST_YEAR,
  LAST_YEAR,
  LAST_EXAM_YEAR,
  MAX_TIME,
  LABS_EXAMS = NULL,
  LABS_OBS = NULL,
  LABS_GRADES = NULL,
  LABS_COV = NULL,
  SUBSET = NULL,
  VERBOSE = FALSE
) {
  GRADES <- check_grades(GRADES)
  colnames(GRADES) <- LABS_EXAMS
  rownames(GRADES) <- LABS_OBS

  TIMES <- check_times(TIMES)
  colnames(TIMES) <- LABS_EXAMS
  rownames(TIMES) <- LABS_OBS

  TODO <- check_todo(TODO)
  colnames(TODO) <- LABS_EXAMS
  rownames(TODO) <- LABS_OBS

  OUTCOME <- check_outcome(OUTCOME)
  names(OUTCOME) <- LABS_OBS

  FIRST_YEAR <- check_year(FIRST_YEAR)
  names(FIRST_YEAR) <- LABS_OBS

  LAST_YEAR <- check_year(LAST_YEAR)
  names(LAST_YEAR) <- LABS_OBS

  LAST_EXAM_YEAR <- check_yearle(LAST_EXAM_YEAR)
  names(LAST_EXAM_YEAR) <- LABS_OBS

  MAX_TIME <- check_maxt(MAX_TIME)
  names(MAX_TIME) <- LABS_OBS

  X <- check_covariates(X)
  colnames(X) <- LABS_COV
  rownames(X) <- LABS_OBS

  n_obs <- nrow(GRADES)
  n_exams <- ncol(GRADES)
  n_cov <- ncol(X)
  n_grades <- length(table(GRADES[, 1]))
  yb <- max(LAST_EXAM_YEAR[LAST_EXAM_YEAR != 100])
  stopifnot(nrow(TIMES) == n_obs)
  stopifnot(ncol(TIMES) == n_exams)
  stopifnot(length(OUTCOME) == n_obs)
  stopifnot(length(FIRST_YEAR) == n_obs)
  stopifnot(length(LAST_YEAR) == n_obs)
  stopifnot(length(LAST_EXAM_YEAR) == n_obs)
  stopifnot(nrow(X) == n_obs)

  if (is.null(SUBSET)) {
    SUBSET <- 1:n_obs
  }

  # parameter dimensions
  dim_grtcm <- n_exams * (n_grades + 3)
  dim_lat <- 2 + 2 * n_cov
  dim_cr <- (yb + 2) * 2 + 1

  # joint checks
  if (!all(which(is.na(GRADES)) == which(is.na(TIMES)))) {
    stop(
      "For one or more exam, some observations have observed grade but missing time or viceversa."
    )
  }

  if (!all(which(!is.na(GRADES)) %in% which(TODO))) {
    stop(
      "Some observations have observed grade on one or more exams not included in the study plan."
    )
  }

  tmp <- sapply(1:n_obs, function(IDX) {
    if (!all(!is.na(GRADES[IDX, TODO[IDX, ]])) & OUTCOME[IDX] == 1) {
      if (VERBOSE) {
        warning(paste0(
          "Student at row ",
          IDX,
          " (ID:",
          LABS_OBS[IDX],
          ") is coded as graduated but has not passed all the exams in his study plan."
        ))
      }
    }
  })
  tmp <- sapply(1:n_obs, function(IDX) {
    if (
      all(!is.na(GRADES[IDX, TODO[IDX, ]])) &
        (OUTCOME[IDX] != 1 & OUTCOME[IDX] != 0)
    ) {
      if (VERBOSE) {
        warning(paste0(
          "Student at row ",
          IDX,
          "(ID:",
          LABS_OBS[IDX],
          ") passed all the exams in his study plan but has outcome ",
          OUTCOME[IDX],
          "."
        ))
      }
    }
  })

  return(list(
    "data_dims" = list(
      n_obs = n_obs,
      n_exams = n_exams,
      n_grades = n_grades,
      n_cov = n_cov,
      yb = yb
    ),
    "par_dims" = list(grtcm = dim_grtcm, lat = dim_lat, cr = dim_cr),
    "labs" = list(
      obs = LABS_OBS,
      exams = LABS_EXAMS,
      grades = LABS_GRADES,
      cov = LABS_COV
    ),
    "todoMat" = TODO[SUBSET, ],
    "gradesMat" = GRADES[SUBSET, ],
    "timeMat" = TIMES[SUBSET, ],
    "outcome" = OUTCOME[SUBSET],
    "first_year" = FIRST_YEAR[SUBSET],
    "last_year" = LAST_YEAR[SUBSET],
    "yle" = LAST_EXAM_YEAR[SUBSET],
    "max_time" = MAX_TIME[SUBSET],
    "X" = as.matrix(X[SUBSET, ], length(SUBSET), n_cov)
  ))
}
