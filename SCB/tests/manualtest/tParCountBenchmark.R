rm (list = ls())

myPath <- doPath()

myTParCountArray  <- 2:40

mySuperReplicationCount <- 20
myReplicationCount <- 20
mySampleSize <- 60
myLag <- 1
myLagCount <- 4
myKernel <- normalDifferenceKernel
myBandwidth <- 0.5
myNonCoverageProbability <- 0.05


fileName <- paste( "Benchmark_TParCount","l", myLag, "bandW", myBandwidth,
                   "sampleS", mySampleSize,
                   "alpha", myNonCoverageProbability, sep = "_")

fileName <- gsub("\\.","", fileName)
myFileName <- gsub(" ","_",paste(fileName ,Sys.time(),".jpg",sep = ""))
myFileName <- gsub(":","_", myFileName)
jpeg(paste(myPath,"/",myFileName,sep=""))


duration <- numeric(length = length(myTParCountArray))

for(i in 1: length(myTParCountArray))
{
  myTParArray  <- createTParArray(tParCount = myTParCountArray[i])

  Start=Sys.time()
  cat("TParCount = ",myTParCountArray[i],"\n")
  nonCoverageFreqArray = computeNonCoverageFreqArray(
    superReplicationCount = mySuperReplicationCount,
    replicationCount = myReplicationCount,
    sampleSize = mySampleSize,
    lag = myLag,
    lagCount = myLagCount,
    tParArray = myTParArray,
    kernel = myKernel,
    bandwidth = myBandwidth,
    nonCoverageProbability = myNonCoverageProbability)

  End=Sys.time()
  duration[i]=difftime(End,Start,units = "secs")
}

plot(x=myTParCountArray,y=duration,main = "Benchmark for TParCount")

myDF <- data.frame(myTParCountArray,duration)

myFileName <- gsub(" ","_",paste(fileName, Sys.time(),".csv",sep = ""))
myFileName <- gsub(":","_",myFileName)
write.csv(myDF, paste(myPath,"/",myFileName,sep=""))

dev.off()

