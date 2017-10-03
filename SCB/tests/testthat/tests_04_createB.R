createBFunction<-function(sampleSize=5)
{
  cat("\n Testing \'tests_04_createB\' \n")
  mySampleSize = sampleSize
  Start=Sys.time()
  b <- computeB(mySampleSize)
  End=Sys.time()
  Duration=End-Start

  cat("bandwidth =", b, "\n")
  cat("Duration =",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("sampleSize= ",sampleSize,"\n")
}
  test_that("Testing  \'createB\'", {
    createBFunction()
    cat(" End of test createB","\n")
    cat("=====================")
  }
  )


