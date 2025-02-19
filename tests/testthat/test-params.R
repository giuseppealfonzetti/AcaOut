n_grades <- 4L
n_exams  <- 3L
n_cov    <- 3L
yb       <- 5L

dim_irt <- n_exams*(n_grades+3)
dim_lat <- 2+2*n_cov
dim_cr  <- 2*(yb+2)+1

labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
labs_cov <- if(n_cov>0) paste0("X",1:n_cov)

set.seed(123)
theta_irt <- rnorm(dim_irt)
theta_lat <- rnorm(dim_lat); theta_lat[2] <- abs(theta_lat[2])
theta_cr <- rnorm(dim_cr)
theta <- c(theta_irt, theta_lat, theta_cr)
mat <- irtVec2Mat(
  THETA_IRT = theta_irt,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades
)

test_that("irtVec2Mat() exams labs",{
  expect_equal(
    rownames(mat),
    labs_exams
  )
})
test_that("irtVec2Mat() grades labs",{
  expect_equal(
    colnames(mat)[1:n_grades],
    labs_grades
  )
})
for (exam in 1:nrow(mat)) {

  test_that("irtVec2Mat() intercepts",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 2, exam-1),
      as.numeric(mat[exam, 1:n_grades])
    )
  })

  test_that("irtVec2Mat() slope",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 1, exam-1),
      as.numeric(mat[exam, n_grades+1])
    )
  })

  test_that("irtVec2Mat() time loc",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 3, exam-1),
      as.numeric(mat[exam, n_grades+2])
    )})

  test_that("irtVec2Mat() time scale",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 4, exam-1),
      as.numeric(mat[exam, n_grades+3])
    )})
}
test_that("irtMat2Vec()", {
  expect_equal(theta_irt, irtMat2Vec(mat))
})
theta[23]
par_list <- parVec2List(
  THETA=theta,
  N_GRADES=n_grades, N_EXAMS=n_exams, N_COV=n_cov, YB=yb,
  LABS_EXAMS=labs_exams, LABS_GRADES=labs_grades, LABS_COV=labs_cov)

test_that("parVec2List/List2Vec",{
  expect_equal(
    parList2Vec(par_list),
    theta
  )
})


test_that("vector reparameterisation",{
  expect_equal(
    parVec2Repar(
      THETA=theta,
      N_GRADES=n_grades, N_EXAMS=n_exams, N_COV=n_cov, YB=yb,
      LABS_EXAMS=labs_exams, LABS_GRADES=labs_grades, LABS_COV=labs_cov, TIDY = TRUE)$par,
    parVec2Repar(
      THETA=theta,
      N_GRADES=n_grades, N_EXAMS=n_exams, N_COV=n_cov, YB=yb,
      LABS_EXAMS=labs_exams, LABS_GRADES=labs_grades, LABS_COV=labs_cov, TIDY = FALSE)
  )
})


