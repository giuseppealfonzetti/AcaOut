#' Fit conditional competing risk model
#'
#' @param DATA data as per [check_data()].
#' @param PAR_START Starting parameter vector
#' @param LATMAT Latent scores to be used for conditional competing risk model
#'
#' @export
fit_CCR <- function(DATA, PAR_START, LATMAT) {
  dim_irt <- DATA$par_dims$grtcm
  dim_lat <- DATA$par_dims$lat

  CCRfun <- function(PAR, GRFLAG) {
    obj <- CCR(
      THETA = c(PAR_START[c(1:(dim_irt + dim_lat))], PAR),
      OUTCOME = DATA$outcome,
      COVARIATES = DATA$X,
      YEAR_FIRST = DATA$first_year,
      YEAR_LAST = DATA$last_year,
      YEAR_LAST_EXAM = DATA$yle,
      LATMAT = as.matrix(LATMAT),
      YB = DATA$data_dims$yb,
      GRFLAG = GRFLAG
    )

    if (GRFLAG) {
      return(-obj$gr[-c(1:(DATA$par_dims$grtcm + DATA$par_dims$lat))])
    } else {
      return(-obj$ll)
    }
  }

  CCRNLL <- function(PAR) {
    CCRfun(PAR, GRFLAG = FALSE)
  }
  CCRNGR <- function(PAR) {
    CCRfun(PAR, GRFLAG = TRUE)
  }

  fit <- ucminf::ucminf(
    par = PAR_START[-c(1:(dim_irt + dim_lat))],
    fn = CCRNLL,
    gr = CCRNGR,
    hessian = 2
  )

  return(
    list(
      "data" = DATA,
      "mod" = "ccr",
      "fit" = fit
    )
  )
  return(fit)
}
