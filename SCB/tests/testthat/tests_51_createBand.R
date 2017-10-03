createBandFunction <- function()
{
  cat("\n Testing \'tests_51_createBand\' \n ")
  mySampleSize = 10
  myTParCount = 10
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockTVMA1Array <-
    createTVMA1CoefArray(sampleSize = mySampleSize)
  myLag = 5
  myLagCount = computeLagCount(lag = myLag, sampleSize = mySampleSize)
  myKernel = normalDifferenceKernel
  myNonCoverageProbability = 0.05
  errorIfNotInputCompatible(alpha = myNonCoverageProbability,lag = myLag,sampleSize = mySampleSize)
  Start = Sys.time()
  mySample <- createSample(sampleSize = mySampleSize)
  band = createBand(
    tParArray = mockTParArray,
    lag = myLag,
    lagCount = myLagCount,
    sampleSize = mySampleSize,
    nonCoverageProbability = myNonCoverageProbability
  )
  End = Sys.time()
  Duration = End - Start


  cat("\nband1", band)
  cat("\nDuration= ", Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("SampleSize= ",mySampleSize,"\n")
  cat("TParCount= ",myTParCount,"\n")
  cat("Lag= ",myLag,"\n")
  cat("LagCount= ",myLagCount,"\n")
  cat("NonCoverageProbability= ",myNonCoverageProbability,"\n")

}

test_that("Testing \'createBand\'", {
  createBandFunction()
  cat("End of test createBand", "\n")
  cat("=====================")
})
