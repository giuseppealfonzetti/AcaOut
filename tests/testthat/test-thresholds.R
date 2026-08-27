test_that("reparThr() output",{
  nthr <- 5

  # reference implementation
  mic_fun <-function(param)
  {
    ncat<-length(param)
    tmp<-param[-ncat]
    tmp[-1]<-exp(tmp[-1])
    param[-ncat]<-cumsum(tmp)
    param
  }

  set.seed(123)
  for (trial in 1:3) {
    param <- rnorm(nthr)
    thr <- reparThr(param, CON2UN = F)

    expect_true(sum(sort(thr, decreasing = F) == thr)==length(thr))

    expect_equal(thr, mic_fun(c(param, 1))[1:nthr])

    expect_equal(reparThr(thr, CON2UN = T), param)

  }
})
