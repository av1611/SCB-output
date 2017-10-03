rm (list = ls())

myPath <- doPath()

myTParCount  <-2
myTParArray  <- createTParArray(tParCount = myTParCount)
mySuperReplicationCount <- 4
myReplicationCount <- 4
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


sampleSizeArray <- c(10,20,30,40,50,60)

duration <- numeric(length = length(sampleSizeArray))

for(i in 1: length(sampleSizeArray))
{
  Start=Sys.time()
  cat("Sample Size = ",sampleSizeArray[i],"\n")
  nonCoverageFreqArray = computeNonCoverageFreqArray(
    superReplicationCount = mySuperReplicationCount,
    replicationCount = myReplicationCount,
    sampleSize = sampleSizeArray[i],
    lag = myLag,
    lagCount = myLagCount,
    tParArray = myTParArray,
    kernel = myKernel,
    bandwidth = myBandwidth,
    nonCoverageProbability = myNonCoverageProbability)

  End=Sys.time()
  duration[i]=End-Start
}

plot(x=sampleSizeArray,y=duration,main = "Benchmark for sampleSize")

myDF <- data.frame(sampleSizeArray,duration)

myFileName <- gsub(" ","_",paste(fileName, Sys.time(),".csv",sep = ""))
myFileName <- gsub(":","_",myFileName)
write.csv(myDF, paste(myPath,"/",myFileName,sep=""))

dev.off()


