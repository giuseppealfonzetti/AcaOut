#' Student careers
#'
#' Dataset collecting information about students careers, including exams-related grades and times as well as academic outcomes.
#'
#' @format ## `stud_careers``
#' \describe{
#'   \item{data_dims}{Gather info about data dimensions. List containing number of students (`n_obs`), number of exams (`n_exams`), number of grades categories (`n_grades`), number of covariates (`n_cov`).}
#'   \item{par_dims}{Gather info about model dimension. List containing number of parameters releted to measurement model (`grtcm`), number of parameters related to outcome model (`cr`), number of parameters related to structural model (`lat`).}
#'   \item{labs}{List containing students, exams and grades labels.}
#'   \item{todoMat}{`n_obs` x `n_exams` matrix of logical values denoting the plan of study for each student.}
#' \item{gradesMat}{`n_obs` x `n_exams` matrix collecting students' grades (categories) on each exam.}
#' \item{timeMat}{`n_obs` x `n_exams` matrix of times (0 denotes enrollment day) before passing each exam.}
#' \item{outcome}{`n_obs`-dimensional vector collecting academic outcome information. Encoding: 0 for still enrolled; 1 for graduation, 2 for dropout; 3 for transfer.}
#' \item{first_year}{`n_obs`-dimensional vector collecting first enrollment year.}
#' \item{last_year}{`n_obs`-dimensional vector collecting last enrollment year.}
#' \item{yle}{`n_obs`-dimensional vector collecting the year of study plan completio.}
#' \item{max_time}{`n_obs`-dimensional vector collecting max observed day.}
#' \item{X}{`n_obs` x `n_exams` matrix of covariates.}
#' }
#'
#' @source Alfonzetti, G., and Battauz, M. "A joint model for academic outcomes and students' exam performance with informative censoring".
"stud_careers"
