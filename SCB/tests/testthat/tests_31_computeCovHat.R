computeCovHatFunction = function ( ) {
  cat("\n Testing \'tests_31_computeCovHat\' \n")
  tParCount=10
  tParArray=createTParArray(tParCount = tParCount)
  mySampleSize=10
  mockTVMA1CoefArray <- createTVMA1CoefArray(sampleSize = mySampleSize)
  mockSample <- createSample(sampleSize = mySampleSize)


  myLag = 1
  errorIfNotInputCompatible(alpha = 0.1,sampleSize = mySampleSize,lag = myLag)
  myKernel = normalDifferenceKernel

  Start=Sys.time()
  mockCovHat <- computeCovHat(tParArray = tParArray,
                              lag = myLag,
                              sample = mockSample)
  End=Sys.time()
  Duration=End-Start

  cat ("mockCovHat[1:5] = ", mockCovHat[1:5], "\n")
  cat ("length(mockCovHat) = ", length(mockCovHat), "\n")
  cat("Duration= ",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("sampleSize= ",mySampleSize,"\n")
  cat("tParCount= ",tParCount,"\n")

  expect_that(mockCovHat, is_a("numeric"))
  # expect_that(dim(mockBand)[1], equals(2))  # the number of rows
  # expect_that(dim(mockBand)[2], equals(10)) # the number of cols

}

test_that("Testing \'computeCovHat\'", {
  computeCovHatFunction()
  cat(" End of test computeCovHat","\n")
  cat("=====================")

}

)
