createTVMA1CoefArrayFunction <- function(samplesize=5) {
  cat("\n Testing \'tests_02_createTVMA1CoefArray\' for length and being numeric \n")


  sampleSize = samplesize
  Start=Sys.time()
  mockTVMA1CoefArray <- createTVMA1CoefArray(sampleSize = sampleSize)
  End=Sys.time()
  Duration=End-Start
  cat ("mockTVMA1CoefArray =", mockTVMA1CoefArray[1:sampleSize], "\n",
       "length(mockTVMA1CoefArray) =", length(mockTVMA1CoefArray), "\n")
  cat("Duration= ",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("sampleSize= ",sampleSize,"\n")
 #expect_that(mockTVMA1CoefArray, is_a("numeric"))

}



test_that("Testing tests_createTVMA1CoefArray for length and being numeric", {
  createTVMA1CoefArrayFunction()

  cat(" End of test createTVMA1CoefArray","\n")
  cat("=====================")
  }
)
