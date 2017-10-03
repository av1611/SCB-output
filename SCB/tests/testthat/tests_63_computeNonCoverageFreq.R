  computeNonCoverageFreqFunction <- function()

  {
    cat ("\n Testing \'tests_63_computeNonCoverageFreq\' \n")

    myTParCount  <- 10
    mockTParArray  <- createTParArray(tParCount = myTParCount)
    myReplicationCount <- 5
    mySampleSize <-100
    mySuperReplicationCount <- 10
    myLag <- 1
    myLagCount <- computeLagCount(lag = myLag,sampleSize = mySampleSize)
    trueCorArray <- computeCor(lag = myLag,tParArray = mockTParArray)

    myKernel <- normalDifferenceKernel
    myNonCoverageProbability <- 0.05

    errorIfNotInputCompatible(alpha = myNonCoverageProbability,
                              lag = myLag,sampleSize = mySampleSize)

    Start=Sys.time()

   nonCoverageFreq <- computeNonCoverageFreq(replicationCount = myReplicationCount,
                                            sampleSize = mySampleSize,
                                            lagCount = myLagCount,
                                            lag = myLag,
                                            tParArray = mockTParArray,
                                            corArray = trueCorArray,
                                            nonCoverageProbability = myNonCoverageProbability,
                                            superReplicationCount = mySuperReplicationCount,
                                            fileName = "tests_63_computeNonCoverageFreq")
End=Sys.time()
Duration=End-Start

   cat("NonCoverageFreq:", nonCoverageFreq)
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
     cat("\nEnd of test of computeNonCoverageFreq","\n")
     cat("=====================")
  }

  test_that("computeNonCoverageFreq", {
    computeNonCoverageFreqFunction()
  })

