test_that("simulate_crirt_data() returns a valid data bundle", {
  sim <- simulate_crirt_data(
    N_STUDENTS = 25,
    N_EXAMS = 6,
    N_GRADES = 4,
    MAX_YEAR = 4,
    N_COV = 3,
    SEED = 42
  )

  expect_type(sim, "list")
  expect_true(all(
    c(
      "todoMat",
      "gradesMat",
      "timeMat",
      "outcome",
      "first_year",
      "last_year",
      "yle",
      "max_time",
      "X",
      "params",
      "latent"
    ) %in%
      names(sim)
  ))

  expect_equal(dim(sim$gradesMat), c(25, 6))
  expect_equal(dim(sim$timeMat), c(25, 6))
  expect_equal(dim(sim$todoMat), c(25, 6))
  expect_equal(length(sim$outcome), 25)
  expect_equal(ncol(sim$X), 3)
  expect_equal(sim$data_dims$n_obs, 25)
  expect_equal(sim$yb, 4)

  expect_true(is.matrix(sim$latent))
  expect_equal(dim(sim$latent), c(25, 2))

  expect_true(all(which(is.na(sim$timeMat)) %in% which(is.na(sim$gradesMat))))
  expect_true(all(sim$max_time <= 4 * 365))
})

test_that("simulate_crirt_data()", {
  base_sim <- simulate_crirt_data(
    N_STUDENTS = 100,
    N_EXAMS = 3,
    N_GRADES = 3,
    MAX_YEAR = 3,
    N_COV = 1,
    SEED = 99
  )

  custom_lat <- matrix(
    c(seq(-1, 0.8, length.out = 100), seq(0.5, -0.5, length.out = 100)),
    ncol = 2
  )

  sim_custom <- simulate_crirt_data(
    N_STUDENTS = 100,
    N_EXAMS = 3,
    N_GRADES = 3,
    MAX_YEAR = 3,
    N_COV = 1,
    SEED = 456,
    PARAMS = base_sim$params,
    LATMAT = custom_lat
  )

  expect_equal(sim_custom$params, base_sim$params)
  expect_equal(sim_custom$latent, custom_lat, ignore_attr = TRUE)
  expect_equal(dim(sim_custom$gradesMat), c(100, 3))
  expect_equal(sim_custom$data_dims$n_obs, 100)
  expect_equal(sim_custom$yb, 3)
})
