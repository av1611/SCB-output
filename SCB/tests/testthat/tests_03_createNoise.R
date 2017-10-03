createNoiseFunction <- function (sampleSize=5) {
  cat("\n Testing \'tests_03_createNoise\' for length and being numeric \n")
  set.seed(100)
  mySampleSize = sampleSize
  Start=Sys.time()
  mockNoise <- createNoise(sampleSize = mySampleSize,
                           mean = 0,
                              sd = 1)
  End=Sys.time()
  Duration=End-Start

  set.seed(100)
  mockRnorm <- rnorm(n = mySampleSize,
                     mean = 0,
                     sd = 1)



  cat ("mockNoise[1:5] =", mockRnorm[1:5], "\n")
  cat ("length(mockNoise) =", length(mockNoise), "\n")
  cat("Duration =",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("sampleSize= ",sampleSize,"\n")

  expect_that(mockRnorm, is_a("numeric"))
}

test_that("Testing \'createNoise\' for length and being numeric\'", {
  mockRnorm <- createNoiseFunction()
  cat(" End of test createNoise ","\n")
  cat("=====================")

}
)




