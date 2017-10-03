computeIsCoveredArrayFunction <- function()
{
  cat("\n Testing \'tests_62_computeIsCoveredArray\'\n")

  mySampleSize=10
  myTParCount = 10
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockTVMA1Array <- createTVMA1CoefArray(sampleSize = mySampleSize)
  myReplicationCount=5
  myLag = 1

  myLagCount = computeLagCount(lag = myLag,sampleSize = mySampleSize)
  mySuperReplicationCount <- 3
  myNonCoverageProbability <- 0.05
  errorIfNotInputCompatible(alpha = myNonCoverageProbability,
                            lag = myLag,sampleSize = mySampleSize)

  bandsBrick=createBandsBrick(tParArray=mockTParArray,
                              lag=myLag,
                              lagCount=myLagCount,

                              sampleSize = mySampleSize,
                              nonCoverageProbability = myNonCoverageProbability,
                              replicationCount=myReplicationCount)
  mockCorArray <- computeCor(lag = myLag,
                             tParArray = mockTParArray)


Start=Sys.time()
      isCoveredArray <-
        computeIsCoveredArray(bandsBrick = bandsBrick,
                              corArray = mockCorArray,
                              sampleSize = mySampleSize,
                              lag = myLag,
                              tParCount = myTParCount,
                              replicationCount = myReplicationCount,
                              superReplicationCount = mySuperReplicationCount,
                              nonCoverageProbability = myNonCoverageProbability,
                              fileName = "tests_62_computeIsCoveredArray")
End=Sys.time()
Duration=End-Start

      # expect(5 == 5)
      cat("\nisCoveredArray: ", isCoveredArray)
      cat("\nDuration= ",Duration,"\n")
      cat("\nTest parameters :","\n")
      cat("SampleSize= ",mySampleSize,"\n")
      cat("TParCount= ",myTParCount,"\n")
      cat("Lag= ",myLag,"\n")
      cat("LagCount= ",myLagCount,"\n")
      cat("ReplicationCount= ",myReplicationCount,"\n")
      cat("SuperReplicationCount= ",mySuperReplicationCount,"\n")
      cat("NonCoverageProbability= ",myNonCoverageProbability,"\n")




}

test_that("Testing \'isCoveredArray\'", {
  computeIsCoveredArrayFunction()
  cat("End of test isCoveredArray","\n")
  cat("=====================")
})

