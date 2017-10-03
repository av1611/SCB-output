computeIsCoveredFunction = function()
{
  cat ("\n Testing \'tests_61_computeIsCovered\'\n")

  mySampleSize = 100
  myTParCount = 10
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockTVMA1Array <-
    createTVMA1CoefArray(sampleSize = mySampleSize)
  myReplicationCount=5
  mySuperReplicationCount <- 3
  myCoefFunction=sin

  myLag = 2

  myLagCount =computeLagCount(lag = myLag,sampleSize = mySampleSize)

  myNonCoverageProbability = 0.05
  errorIfNotInputCompatible(alpha = myNonCoverageProbability,lag = myLag,sampleSize = mySampleSize)
  fileName <- paste("ss", mySampleSize, "l", myLag, "alpha"
                    , myNonCoverageProbability, sep = "_")
  mockBand = createBand(
    tParArray = mockTParArray,
    lag = myLag,
    lagCount = myLagCount,

    sampleSize = mySampleSize,
    nonCoverageProbability = myNonCoverageProbability
  )

  mockTVMA1CoefArray <- createTVMA1CoefArray(sampleSize = mySampleSize)

  mockCor <- computeCor(lag = myLag,
                        tParArray = mockTParArray)

  Start=Sys.time()
  isCovered <- computeIsCovered(mockBand,
                                mockCor,
                                sampleSize = mySampleSize,
                                lag = myLag,
                                replicationCount = myReplicationCount,
                                superReplicationCount = mySuperReplicationCount,
                                nonCoverageProbability = myNonCoverageProbability,
                                fileName = "test_61_computeIsCovered")
  End=Sys.time()
  Duration=End-Start


  cat("\nIsCovered:",isCovered, "\n")
  cat("Duration",Duration,"\n")
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

test_that("Testing computeIsCovered", {
  computeIsCoveredFunction()
  cat("End of test computeIsCovered","\n")
  cat("=====================")
})
