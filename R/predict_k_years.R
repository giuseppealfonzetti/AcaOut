#' K-year-ahead outcome predictions
#'
#' Compute the probability that each student experiences graduation, dropout,
#' transfer, or remains enrolled by the end of academic year
#' `CUTOFF_YEAR + HORIZON`, conditioning on the data available up to the end of
#' `CUTOFF_YEAR`. Probabilities are obtained by iterating
#' [predict_one_year()] over successive horizons and accumulating the implied
#' hazards.
#'
#' @param FIT Fitted object returned by [fit_EM()] or [fit_BFGS()] with
#'   `mod = "full"`.
#' @param CUTOFF_YEAR Integer year `y` whose information is used as the starting
#'   point. Must satisfy `0 <= CUTOFF_YEAR`.
#' @param HORIZON Strictly positive integer defining how many future academic
#'   years are evaluated. The final horizon year is `CUTOFF_YEAR + HORIZON`.
#' @inheritParams predict_one_year
#'
#' @return A tibble with one row per student containing cumulative probabilities
#'   of graduation, dropout, transfer, and remaining enrolled by the end of the
#'   target horizon.
#'
#' @export
#'
#' @importFrom tibble tibble
predict_k_years <- function(
  FIT,
  CUTOFF_YEAR,
  HORIZON,
  YEAR_LENGTH = 365,
  GRAD_EXTENSION = 1.5
) {
  if (!(FIT$mod %in% "full")) {
    stop("Prediction is only available for `mod = \"full\"` fits.")
  }
  if (length(CUTOFF_YEAR) != 1 || !is.finite(CUTOFF_YEAR)) {
    stop("`CUTOFF_YEAR` must be a single finite value.")
  }
  CUTOFF_YEAR <- as.integer(CUTOFF_YEAR)
  if (CUTOFF_YEAR < 0) {
    stop("`CUTOFF_YEAR` must be non-negative.")
  }
  if (length(HORIZON) != 1 || !is.finite(HORIZON)) {
    stop("`HORIZON` must be a single finite value.")
  }
  HORIZON <- as.integer(HORIZON)
  if (HORIZON <= 0) {
    stop("`HORIZON` must be strictly positive.")
  }
  if (length(YEAR_LENGTH) != 1 || !is.finite(YEAR_LENGTH) || YEAR_LENGTH <= 0) {
    stop("`YEAR_LENGTH` must be a positive scalar.")
  }
  if (length(GRAD_EXTENSION) != 1 || !is.finite(GRAD_EXTENSION) || GRAD_EXTENSION <= 1) {
    stop("`GRAD_EXTENSION` must be greater than 1.")
  }

  yb <- FIT$data$data_dims$yb
  if (CUTOFF_YEAR + HORIZON > yb) {
    stop("`CUTOFF_YEAR + HORIZON` exceeds the maximum modelled year.")
  }

  subject_id <- FIT$data$labs$obs
  n_students <- length(subject_id)
  survival <- rep(1, n_students)
  cum_grad <- cum_drop <- cum_trans <- rep(0, n_students)

  current_cutoff <- CUTOFF_YEAR
  for (step in seq_len(HORIZON)) {
    res <- suppressMessages(
      predict_one_year(
        FIT,
        CUTOFF_YEAR = current_cutoff,
        YEAR_LENGTH = YEAR_LENGTH,
        GRAD_EXTENSION = GRAD_EXTENSION
      )
    )

    res <- res[match(subject_id, res$subject_id), ]
    if (anyNA(res$prob_graduation) || anyNA(res$prob_dropout) ||
        anyNA(res$prob_transfer) || anyNA(res$prob_still_enrolled)) {
      stop("Unexpected NA probabilities returned by `predict_one_year()`.")
    }

    g <- res$prob_graduation
    d <- res$prob_dropout
    tr <- res$prob_transfer
    total <- pmin(pmax(g + d + tr, 0), 1)

    cum_grad <- cum_grad + survival * g
    cum_drop <- cum_drop + survival * d
    cum_trans <- cum_trans + survival * tr

    survival <- survival * (1 - total)
    survival <- pmin(pmax(survival, 0), 1)

    current_cutoff <- current_cutoff + 1L
  }

  prob_still <- survival

  occurred_before_cutoff <- FIT$data$outcome %in% c(1L, 2L, 3L) &
    FIT$data$last_year <= CUTOFF_YEAR
  if (any(occurred_before_cutoff)) {
    grad_idx <- occurred_before_cutoff & FIT$data$outcome == 1L
    drop_idx <- occurred_before_cutoff & FIT$data$outcome == 2L
    trans_idx <- occurred_before_cutoff & FIT$data$outcome == 3L
    cum_grad[grad_idx] <- 1
    cum_drop[drop_idx] <- 1
    cum_trans[trans_idx] <- 1
    prob_still[occurred_before_cutoff] <- 0
  }

  tibble(
    subject_id = subject_id,
    cutoff_year = CUTOFF_YEAR,
    horizon = HORIZON,
    target_year = CUTOFF_YEAR + HORIZON,
    prob_graduation = pmin(pmax(cum_grad, 0), 1),
    prob_dropout = pmin(pmax(cum_drop, 0), 1),
    prob_transfer = pmin(pmax(cum_trans, 0), 1),
    prob_still_enrolled = pmin(pmax(prob_still, 0), 1)
  )
}
