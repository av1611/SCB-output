saveDoubleAplhaHatArrayFunction <- function()
{
  cat ("\n Testing \'tests_73_saveDoubleAplhaHatArray\'\n")
  myTParCount  <- 10
  myTParArray  <- createTParArray(tParCount = myTParCount)
  mySuperReplicationCount <- 1
  myReplicationCount <- 42
  sampleSize <- 5
  myLag <- 2
  myLagCount <- computeLagCount(lag = myLag,sampleSize = sampleSize)
  nonCoverageProbabilities <- c(0.2,0.4,0.6,0.8)

  alphaHats<-createDoubleAlphaArray(superReplicationCount = mySuperReplicationCount,
                         replicationCount = myReplicationCount,
                         sampleSize = sampleSize,
                         lag = myLag, lagCount = myLagCount,
                         alphaArray = nonCoverageProbabilities,
                         tParArray = myTParArray,
                         fileName = "tests_73_SaveDoubleAlphaHatArray")

   saveDoubleAplhaHatArray(nonCoverageProbabilities = nonCoverageProbabilities,
                          alphaHats = alphaHats, sampleSize = sampleSize,
                          lag = myLag, replicationCount = myReplicationCount,
                          superReplicationCount = mySuperReplicationCount,
                          fileName = "tests_73_SaveDoubleAlphaHatArray")
   cat("\n=====================")
   cat("\nTest parameters :","\n")
   cat("SampleSize= ",sampleSize,"\n")
   cat("TParCount= ",myTParCount,"\n")
   cat("Lag= ",myLag,"\n")
   cat("LagCount= ",myLagCount,"\n")
   cat("NonCoverageProbabilities= ",nonCoverageProbabilities,"\n")
   cat("ReplicationCount= ",myReplicationCount,"\n")
   cat("SuperReplicationCount= ",mySuperReplicationCount,"\n")


 }
test_that("Testing \'saveNonCoverageFreqArray\'", {
  saveDoubleAplhaHatArrayFunction()
  cat("\nEnd of test saveNonCoverageFreqArray", "\n")
  cat("=====================")
})

