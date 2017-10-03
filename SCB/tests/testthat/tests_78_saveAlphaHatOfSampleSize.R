saveAlphaHatOfSampleSizeFunction <- function()
{
  cat ("\n Testing \'tests_78_saveAlphaHatOfSampleSize\'\n")
  sampleSize <- seq(10,100,by=10)

  nonCoverageProbability <- 0.2

  tParCount <- 10
  tParArray <- createTParArray(tParCount = tParCount)

  lag <- 2

  kernel <- normalDifferenceKernel

  replicationCount <- 42
  mySuperReplicationCount <- 1

  alphaHat <- matrix(0,nrow = mySuperReplicationCount,ncol = length(sampleSize))
  for(index in 1:length(sampleSize))
  {
    cat("\nsampleSize = ",sampleSize[index])
    bandwidth <- computeB(sampleSize[index])
    lagCount <- computeLagCount(sampleSize[index],lag)

    errorIfNotInputCompatible(alpha = nonCoverageProbability,
                              lag = lag,
                              sampleSize = sampleSize[index])

    alphaHat[,index] <- computeNonCoverageFreqArray(
                           superReplicationCount = mySuperReplicationCount,
                           replicationCount = replicationCount,
                           sampleSize = sampleSize[index],
                           lag = lag,
                           lagCount = lagCount,
                           tParArray = tParArray,
                           nonCoverageProbability = nonCoverageProbability,
                           fileName = "tests_78_saveAlphaHatOfSampleSize")
  }
  fileName <- "AlphaOfSampleSize"
  fileName <- paste(fileName, "l", lag, sep = "_")
  fileName <- paste(fileName, "repC", replicationCount,"SrepC", mySuperReplicationCount, sep = "_")

  xLab <- "sampleSize"

  subTitle <- paste(xLab,"\nlag= ",lag, ", tParCount = " ,tParCount,
                    ",\n replicationCount= ", replicationCount, ", 'SuperRep = ",
                    mySuperReplicationCount,sep = "")

  path <- doPath()
  saveJpg(fileName,path)
  df <- rbind(sampleSize,alphaHat)
  saveCVS(fileName,path,df)
  alpha <- rep(nonCoverageProbability,times = length(sampleSize))

  plot(alpha~sampleSize,type="c",ylim = c(0,1),col="blue",xlab="",ylab = "alphaHat")

  title(main = "AlphaHat Of SampleSize", sub = subTitle)

  for(index in 1:length(sampleSize))
  {
    x <- rep(sampleSize[index],times = mySuperReplicationCount)
    points(alphaHat[,index]~x,pch=20)
  }
  graphics.off()
}
test_that("testing saveAlphaHatOfSampleSize",
          {
            saveAlphaHatOfSampleSizeFunction()
            cat(" End of test saveAlphaHatOfSampleSize", "\n")
            cat("=========================================")
          })
