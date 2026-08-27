n_grades <- 4
n_exams <- 3

set.seed(123)
theta_irt <- rnorm(n_exams*(n_grades+3))
mat_irt <- matrix(theta_irt, nrow = n_exams, byrow = TRUE)
for (exam in 0:(n_exams-1)) {
  test_that("extract_params_irt() intercepts",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 2, exam),
      reparThr(mat_irt[exam+1, 1:n_grades])
    )
  })
  test_that("extract_params_irt() slopes",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 1, exam),
      mat_irt[exam+1, n_grades+1]
    )
  })

  test_that("extract_params_irt() time location",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 3, exam),
      mat_irt[exam+1, n_grades+2]
    )
  })

  test_that("extract_params_irt() time scale",{
    expect_equal(
      extract_params_irt(theta_irt, n_grades, n_exams, 4, exam),
      exp(mat_irt[exam+1, n_grades+3])
    )
  })
}

t(apply(mat_irt[,1:n_grades], 1, reparThr, CON2UN = FALSE))
reparThr(mat_irt[3,1:n_grades], CON2UN = FALSE)
