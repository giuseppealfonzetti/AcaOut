#' Get Info about the dataset
#'
#' @description Helper function to extract information needed from the administrative dataset
#'
#' @param DATA_IRT Data for the IRT model
#' @param DATA_CR Data for the CR model
#' @param GRADE_COL Label of the column where to find the grades
#' @param EXAM_COL Label of the column where to find exams ids
#' @param STUDID_COL Label of the column where to find students ids.
#' @param TIME_COL Label of the column where to find exams times.
#' @param MAXTIME_COL Label of the column where to find the maximum time observed for each student.
#'
#' @export
#' @importFrom rlang .data
#' @importFrom dplyr select arrange mutate group_by starts_with ungroup across everything left_join rename if_else
#' @importFrom tidyr pivot_wider nest unnest replace_na
dataRestruct <- function(
  DATA_IRT = NULL,
  DATA_CR = NULL,
  GRADE_COL = "xcat",
  EXAM_COL = "itemID",
  STUDID_COL = "subjectID",
  TIME_COL = "time",
  MAXTIME_COL = "maxTime"
) {
  grades_labs <- sort(unique(DATA_IRT[, GRADE_COL]), decreasing = FALSE)
  exams_labs <- sort(unique(DATA_IRT[, EXAM_COL]))
  ids <- sort(unique(DATA_IRT[, STUDID_COL]))

  n <- nrow(DATA_IRT)
  n_exams <- length(exams_labs)
  n_grades <- length(grades_labs)

  exam_grades_tab <- table(
    DATA_IRT[, GRADE_COL],
    DATA_IRT[, EXAM_COL],
    useNA = 'ifany'
  )

  nestedIRT <- DATA_IRT |>
    select(
      'subject_id' = STUDID_COL,
      'max_time' = MAXTIME_COL,
      'exam_id' = EXAM_COL,
      'grade' = GRADE_COL,
      'time' = TIME_COL
    ) |>
    arrange(.data$exam_id) |>
    mutate(todo = TRUE, obsflag = !is.na(.data$exam_id)) |>
    pivot_wider(
      id_cols = c('subject_id', 'max_time'),
      names_from = 'exam_id',
      values_from = c('todo', 'grade', 'time', 'todo')
    ) |>
    group_by(.data$subject_id) |>
    nest(
      todo = starts_with('todo'),
      grades = starts_with('grade'),
      times = starts_with('time')
    ) |>
    arrange(.data$subject_id)

  todoMat <- nestedIRT |>
    unnest('todo') |>
    select(.data$subject_id, starts_with('todo')) |>
    arrange(.data$subject_id) |>
    mutate(across(.cols = everything(), ~ replace_na(., FALSE))) |>
    ungroup(.data$subject_id) |>
    select(-.data$subject_id) |>
    as.matrix()

  gradesMat <- nestedIRT |>
    unnest('grades') |>
    select(.data$subject_id, starts_with('grade')) |>
    arrange(.data$subject_id) |>
    ungroup(.data$subject_id) |>
    select(-.data$subject_id) |>
    as.matrix()

  gradesMat_idx <- matrix(
    sapply(gradesMat, function(x) {
      if (!is.na(x)) {
        which(grades_labs == as.character(x))
      } else {
        x
      }
    }),
    nrow(gradesMat),
    ncol(gradesMat)
  )
  gradesMat <- gradesMat_idx

  timeMat <- nestedIRT |>
    unnest('times') |>
    select(.data$subject_id, starts_with('time')) |>
    arrange(.data$subject_id) |>
    ungroup(.data$subject_id) |>
    select(-.data$subject_id) |>
    as.matrix()

  colnames(todoMat) <- colnames(gradesMat) <- colnames(timeMat) <- exams_labs
  rownames(todoMat) <- rownames(gradesMat) <- rownames(timeMat) <- ids

  cr_main <- data.frame(
    subject_id = DATA_CR$subjectID,
    DATA_CR$X,
    LFLAG = DATA_CR$L,
    outcome = DATA_CR$outcome
  ) |>
    mutate(
      subject_id = factor(.data$subject_id, levels = ids, ordered = TRUE)
    ) |>
    arrange((.data$subject_id)) |>
    nest(year = starts_with("factor")) |>
    mutate(
      year = map(
        .data$year,
        ~ {
          mat <- as.matrix(.x)
          as.numeric(mat %*% matrix(1:ncol(mat), ncol(mat), 1))
        }
      ),
      outcome = .data$outcome - 1
    ) |>
    unnest(.data$year)

  crtib <- cr_main |>
    select(-.data$LFLAG) |>
    group_by(.data$subject_id) |>
    filter(.data$outcome == max(.data$outcome)) |>
    filter(.data$year == max(.data$year)) |>
    rename(last_year = .data$year) |>
    left_join(
      cr_main |>
        select(.data$subject_id, .data$year, .data$LFLAG) |>
        group_by(.data$subject_id) |>
        filter(.data$LFLAG == max(.data$LFLAG)) |>
        filter(.data$year == min(.data$year)) |>
        mutate(yle = if_else(.data$LFLAG == 0, NA, .data$year)) |>
        select(.data$subject_id, .data$yle),
      by = "subject_id"
    ) |>
    left_join(
      cr_main |>
        select(-.data$LFLAG) |>
        group_by(.data$subject_id) |>
        filter(.data$year == min(.data$year)) |>
        rename(first_year = .data$year) |>
        select(.data$subject_id, .data$first_year),
      by = "subject_id"
    )

  yb = max(crtib$yle, na.rm = TRUE)

  outcome_vec <- crtib$outcome
  fy_vec <- crtib$first_year
  ly_vec <- crtib$last_year
  yle_vec <- if_else(is.na(crtib$yle), 100, crtib$yle)
  X <- as.matrix(
    crtib |>
      ungroup() |>
      select(
        -c(
          .data$subject_id,
          .data$outcome,
          .data$last_year,
          .data$first_year,
          .data$yle
        )
      )
  )
  n_cov = ncol(X)
  fulldata <- nestedIRT |> left_join(crtib, by = "subject_id")

  dim_irt <- n_exams * (n_grades + 3)
  dim_lat <- 2 + 2 * n_cov
  dim_cr <- (yb + 2) * 2 + 1

  out <- list(
    "obs_ids" = ids,
    "exams_tab" = exam_grades_tab,
    "grades_labs" = grades_labs,
    "exams_labs" = exams_labs,
    "n_observations" = n,
    "n_exams" = n_exams,
    "n_grades" = n_grades,
    "n_cov" = n_cov,
    "yb" = yb,
    "cr_ext_cov" = n_cov,
    "todoMat" = todoMat,
    "gradesMat" = gradesMat,
    "timeMat" = timeMat,
    "dim_irt" = dim_irt,
    "dim_cr" = dim_cr,
    "dim_lat" = dim_lat,
    "outcome" = outcome_vec,
    "first_year" = fy_vec,
    "last_year" = ly_vec,
    "yle" = yle_vec,
    "X" = X,
    "fulldata" = fulldata
  )

  return(out)
}
