rm (list = ls())

myPath <- doPath()

myTParCount  <-2
myTParArray  <- createTParArray(tParCount = myTParCount)
mySuperReplicationCount <- 20
myReplicationCount <- 20
mySampleSize <- 50
myLag <- 1
myLagCount <- 4
myKernel <- normalDifferenceKernel
myBandwidth <- 0.5
myNonCoverageProbability <- 0.05


fileName <- paste( "Benchmark_Samplesize","l", myLag, "bandW", myBandwidth, "alpha"
                   , myNonCoverageProbability, sep = "_")

fileName <- gsub("\\.","", fileName)
myFileName <- gsub(" ","_",paste(fileName ,Sys.time(),".jpg",sep = ""))
myFileName <- gsub(":","_", myFileName)
jpeg(paste(myPath,"/",myFileName,sep=""))


mySuperReplicationCount <- 4:20

duration <- numeric(length = length(mySuperReplicationCount))

for(i in 1: length(mySuperReplicationCount))
{
  Start=Sys.time()
  cat("Super Replication count = ",mySuperReplicationCount[i],"\n")
  nonCoverageFreqArray = computeNonCoverageFreqArray(
    superReplicationCount = mySuperReplicationCount[i],
    replicationCount = myReplicationCount,
    sampleSize = mySampleSize,
    lag = myLag,
    lagCount = myLagCount,
    tParArray = myTParArray,
    kernel = myKernel,
    bandwidth = myBandwidth,
    nonCoverageProbability = myNonCoverageProbability)

  End=Sys.time()
  duration[i]=End-Start
}

plot(x=mySuperReplicationCount,y=duration,main = "Benchmark for Super Replication Count")

myDF <- data.frame(mySuperReplicationCount,duration)

myFileName <- gsub(" ","_",paste(fileName, Sys.time(),".csv",sep = ""))
myFileName <- gsub(":","_",myFileName)
write.csv(myDF, paste(myPath,"/",myFileName,sep=""))

dev.off()

