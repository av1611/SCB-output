createBandsBrickFunction <- function()
{
  cat ("\n Testing \'tests_52_createBandsBrick\'\n")
  mySampleSize=40
  myTParCount = 2
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockTVMA1Array <- createTVMA1CoefArray(sampleSize = mySampleSize)
  myReplicationCount=5
  myLag = 1
  myLagCount = computeLagCount(lag = myLag,sampleSize = mySampleSize)
  myNonCoverageProbability=0.05
  errorIfNotInputCompatible(alpha = myNonCoverageProbability,lag = myLag,sampleSize = mySampleSize)
Start=Sys.time()
bandsBrick=createBandsBrick(tParArray=mockTParArray,
                            lag=myLag,
                            lagCount=myLagCount,
                            sampleSize=mySampleSize,
                            nonCoverageProbability = myNonCoverageProbability,
                            replicationCount=myReplicationCount)
End=Sys.time()
Duration=End-Start
cat("bandsArray", bandsBrick[1:5])
cat("\nDuration",Duration,"\n")
cat("=====================")
cat("\nTest parameters :","\n")
cat("SampleSize= ",mySampleSize,"\n")
cat("TParCount= ",myTParCount,"\n")
cat("Lag= ",myLag,"\n")
cat("LagCount= ",myLagCount,"\n")
cat("ReplicationCount= ",myReplicationCount,"\n")
cat("NonCoverageProbability= ",myNonCoverageProbability,"\n")
}
test_that("createBandsBrick", {
  createBandsBrickFunction()
  cat("End of test createBandsBrick","\n")
  cat("=====================")
})
