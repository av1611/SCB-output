saveNonCoverageFreqArrayFunction <- function()
{
  cat ("\n Testing \'tests_72_saveNonCoverageFreqArray\'\n")

  myTParCount  <-10
  myTParArray  <- createTParArray(tParCount = myTParCount)
  mySuperReplicationCount <- 7
  myReplicationCount <- 6
  mySampleSize <- 200
  myLag <- 1
  myLagCount <- computeLagCount(lag = myLag,sampleSize = mySampleSize)
  myNonCoverageProbability <- 0.05

  errorIfNotInputCompatible(alpha = myNonCoverageProbability,
                            lag = myLag,sampleSize = mySampleSize)

  isCompatible = isLagCompatibleWithSampleSize (sampleSize = mySampleSize, lag = myLag)
  cat ("\n isCompatible = ", isCompatible)
  if (! isCompatible) cat ("sampleSize is not compatible with lag, stop")
Start=Sys.time()
  nonCoverageFreqArray = computeNonCoverageFreqArray(
    superReplicationCount = mySuperReplicationCount,
    replicationCount = myReplicationCount,
    sampleSize = mySampleSize,
    lag = myLag,
    lagCount = myLagCount,
    tParArray = myTParArray,
    nonCoverageProbability = myNonCoverageProbability,
    fileName = "tests_72_saveNonCoverageFreqArray"
  )

  nonCoverageProbability <- myNonCoverageProbability
  alphaHatArray <- nonCoverageFreqArray

  saveNonCoverageFreqArray(nonCoverageProbability = nonCoverageProbability,
                                   alphaHatArray = alphaHatArray, sampleSize = mySampleSize,
                           replicationCount = myReplicationCount,
                           lag = myLag, superReplicationCount = mySuperReplicationCount,
                           fileName = "tests_72_saveNonCoverageFreqArray")


  End=Sys.time()
  Duration=End-Start
  cat("\nDuration",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("SampleSize= ",mySampleSize,"\n")
  cat("TParCount= ",myTParCount,"\n")
  cat("Lag= ",myLag,"\n")
  cat("LagCount= ",myLagCount,"\n")
  cat("NonCoverageProbability= ",myNonCoverageProbability,"\n")
  cat("ReplicationCount= ",myReplicationCount,"\n")
  cat("SuperReplicationCount= ",mySuperReplicationCount,"\n")

}

test_that("Testing \'saveNonCoverageFreqArray\'", {
  saveNonCoverageFreqArrayFunction()
  cat("\nEnd of test saveNonCoverageFreqArray", "\n")
  cat("=====================")
})
