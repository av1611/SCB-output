computeNonCoverageFreqArrayFunction <- function() {
  cat("\n Testing \'tests_64_computeNonCoverageFreqArray\' \n")

  myTParCount  <-5
  myTParArray  <- createTParArray(tParCount = myTParCount)
  mySuperReplicationCount <- 7
  myReplicationCount <- 6
  mySampleSize <- 13
  myLag <- 5
  myLagCount <- computeLagCount(lag = myLag,sampleSize = mySampleSize)
  myNonCoverageProbability <- 0.05

  errorIfNotInputCompatible(alpha = myNonCoverageProbability,
                            lag = myLag,sampleSize = mySampleSize)
Start=Sys.time()
  nonCoverageFreqArray <- computeNonCoverageFreqArray(
    superReplicationCount = mySuperReplicationCount,
    replicationCount = myReplicationCount,
    sampleSize = mySampleSize,
    lag = myLag,
    lagCount = myLagCount,
    tParArray = myTParArray,
    nonCoverageProbability = myNonCoverageProbability,
    fileName = "tests_64_computeNonCoverageFreqArray"
  )
  End=Sys.time()
  Duration=End-Start
cat("\n NonCoverageFreqArray= ",nonCoverageFreqArray)

  cat("\nDuration= ",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("SampleSize= ",mySampleSize,"\n")
  cat("TParCount= ",myTParCount,"\n")
  cat("Lag= ",myLag,"\n")
  cat("LagCount= ",myLagCount,"\n")
  cat("ReplicationCount= ",myReplicationCount,"\n")
  cat("SuperReplicationCount= ",mySuperReplicationCount,"\n")
  cat("NonCoverageProbability= ",myNonCoverageProbability,"\n")


}

test_that("Testing \'computeNonCoverageFreqArray\'", {
  computeNonCoverageFreqArrayFunction()
  cat("\nEnd of test computeNonCoverageFreqArray", "\n")
  cat("=====================")
})
