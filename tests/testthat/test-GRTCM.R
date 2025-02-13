n <- 10
set.seed(123)

### gen params ####
n_grades <- 4L
n_exams <- 10L
yb <- 5
n_cov <- 1
dim_cr <- 2*(yb+n_cov+2)+1
dim_irt_lat <- n_exams*(n_grades+3)+2
labs_exams <- paste0('ECO0',1:n_exams)
labs_grades <- c('[18,22)', '[22,25)', '[25,28)', '[29,30L]')
theta <- c(rnorm(dim_irt_lat), rep(NA, dim_cr))
parList <- parVec2List(
  THETA = theta,
  N_GRADES = n_grades,
  N_EXAMS = n_exams,
  N_COV = n_cov,
  YB=yb,
  LABS_EXAMS = labs_exams,
  LABS_GRADES = labs_grades)
irtMat <- parList$irt
# theta_irt <- irtMat2Vec(irtMat)

latMat <- matrix(rnorm(n*2), n, 2)
nodes <- expand.grid(c(1,-1),c(1,-1))
weights <- rep(.25, 4)

#### sim grades ####
gradesMat <- matrix(0, n, n_exams)

for (i in 1:n) {
  for (e in 1:n_exams) {
    #linear predictor exams X grades
    linp <- irtMat[e, n_grades+1] * latMat[i,1] - irtMat[e, 1:n_grades]

    # probabilities of greater grades
    pgg <- exp(linp)/(1+exp(linp))

    # probabilities of grades
    pg <- c(pgg[1:(n_grades-1)] - pgg[2:n_grades], pgg[n_grades])
    gradesMat[i,e] <- which(rmultinom(n = 1, size=1, prob = c(1-sum(pg), pg))==1)-1

  }
}



#### sim times ####
set.seed(123)
timeMat <- matrix(0, n, n_exams)
for (i in 1:n) {
  for (e in 1:n_exams) {
    timeMat[i,e] <- exp(
      rnorm(1,
            mean = irtMat[e, n_grades + 2]-latMat[i,2],
            sd = 1/irtMat[e, n_grades + 3])
    )
  }
}
timeMat[gradesMat==0] <- NA

#### mat to-do ####
todoMat <- matrix(1, n, n_exams)

#### censoring ####
max_day <- max(timeMat, na.rm = TRUE)+10
timeMat[timeMat>max_day] <- NA
obsMat <- matrix(1, n, n_exams)
obsMat[is.na(timeMat)] <- 0


##### TESTS ####
RFUN <- function(x, ROTATE){
  GRTCM_GH(
    THETA = x,
    EXAMS_GRADES = gradesMat,
    EXAMS_DAYS = timeMat,
    EXAMS_SET = todoMat,
    EXAMS_OBSFLAG = obsMat,
    MAX_DAY = rep(max_day,n),
    GRID = as.matrix(nodes),
    WEIGHTS = weights,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    GRFLAG = FALSE,
    ROTGRID = ROTATE
  )$ll
}

RFUN(theta[1:dim_irt_lat], TRUE)

#### check gradient
test_that("Check gradient derivative without node rotation", {


  fit <- GRTCM_GH(
    THETA = theta,
    EXAMS_GRADES = gradesMat,
    EXAMS_DAYS = timeMat,
    EXAMS_SET = todoMat,
    EXAMS_OBSFLAG = obsMat,
    MAX_DAY = rep(max_day,n),
    GRID = as.matrix(nodes),
    WEIGHTS = weights,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    GRFLAG = TRUE,
    ROTGRID = FALSE
  )

  expect_equal(
    numDeriv::grad(RFUN, x = theta[1:dim_irt_lat], ROTATE = FALSE),
    fit$gr[1:dim_irt_lat]
  )
})

test_that("Check gradient derivative with node rotation", {

  fit <- GRTCM_GH(
    THETA = theta,
    EXAMS_GRADES = gradesMat,
    EXAMS_DAYS = timeMat,
    EXAMS_SET = todoMat,
    EXAMS_OBSFLAG = obsMat,
    MAX_DAY = rep(max_day,n),
    GRID = as.matrix(nodes),
    WEIGHTS = weights,
    N_GRADES = n_grades,
    N_EXAMS = n_exams,
    GRFLAG = TRUE,
    ROTGRID = TRUE
  )

  expect_equal(
    numDeriv::grad(RFUN, x = theta[1:dim_irt_lat], ROTATE = TRUE),
    fit$gr[1:dim_irt_lat]
  )
})

