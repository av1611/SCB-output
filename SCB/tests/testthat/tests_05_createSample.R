createSampleFunction <- function (sampleSize=5)
{
  cat("\n Testing \'tests_05_createSample\' for length and being numeric \n")
  mySampleSize =sampleSize
  Start=Sys.time()
  mockSample <- createSample(sampleSize = mySampleSize)
  End=Sys.time()
  Duration=End-Start

  cat ("mockSample[1:5] =", mockSample[1:5], "\n")
  cat("Duration= ",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("sampleSize= ",sampleSize,"\n")
  expect_that(mockSample, is_a("array"))
}



test_that("Testing \'createSample\' for length and being numeric",
{
  createSampleFunction()
  cat(" End of test createSample","\n")
  cat("=====================")
}
)
