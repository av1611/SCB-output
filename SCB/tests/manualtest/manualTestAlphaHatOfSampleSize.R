sampleSize <- seq(10,200,by=50)

nonCoverageProbability <- 0.2

tParCount <- 10
tParArray <- createTParArray(tParCount = tParCount)

lag <- 2

kernel <- normalDifferenceKernel

replicationCount <- 6
mySuperReplicationCount <- 7

alphaHat <- matrix(0,nrow = mySuperReplicationCount,ncol = length(sampleSize))
for(index in 1:length(sampleSize))
{
  bandwidth <- computeB(sampleSize[index])
  lagCount <- computeLagCount(sampleSize[index],lag)

  alphaHat[,index] <- computeNonCoverageFreqArray(
                         superReplicationCount = mySuperReplicationCount,
                         replicationCount = replicationCount,
                         sampleSize = sampleSize[index],
                         lag = lag,
                         lagCount = lagCount,
                         tParArray = tParArray,
                         kernel = kernel,
                         bandwidth = bandwidth,
                         nonCoverageProbability = nonCoverageProbability)
}
fileName <- "AlphaOfSampleSize"
fileName <- paste(fileName, "ss", sampleSize, "l", lag, "bandW", bandwidth, sep = "_")
fileName <- paste(fileName, "repC", replicationCount,"SrepC", mySuperReplicationCount, sep = "_")

xLab <- "sampleSize"

subTitle <- paste(xLab,"\nlag= ",lag, ", tParCount = " ,tParCount,
                  ",\n replicationCount= ", replicationCount, ", 'SuperRep = ",
                  mySuperReplicationCount,sep = "")

path <- doPath()
saveJpg(fileName,path)
#df <- data.frame(alphaHat)
#saveCVS(fileName,path,df)
alpha <- rep(nonCoverageProbability,times = length(sampleSize))

plot(alpha~sampleSize,type="c",ylim = c(0,1),col="blue",xlab="",ylab = "alphaHat")

title(main = "Alpha Of SampleSize", sub = subTitle)

for(index in 1:length(sampleSize))
{
  x <- rep(sampleSize[index],times = mySuperReplicationCount)
  points(jitter(alphaHat[,index])~x,pch=20)
}
graphics.off()

