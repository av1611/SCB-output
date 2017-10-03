saveMeForSampleSizeFunction <- function()
{
  cat ("\n Testing \'tests_81_saveMeForSampleSize\'\n")
  sampleSize <- seq(10,150,by=100)
  tParCount <- 10
  tParArray <- createTParArray(tParCount)
  tParArrayCenter <- round(length(tParArray)/2,0)
  nonCoverageProbability <- 0.2
  lag <- 2
  nMe <- 5
  meCenterArray <- matrix(0,nrow = length(sampleSize),ncol = nMe)
  for(indexCol in 1:nMe)
  {
    for(indexRow in 1:length(sampleSize))
    {
      errorIfNotInputCompatible(alpha = nonCoverageProbability,lag = lag,
                                sampleSize = sampleSize[indexRow])
      lagCount <- computeLagCount(sampleSize = sampleSize[indexRow],lag = lag)
      sample <- createSample(sampleSize = sampleSize[indexRow])
      allCor <- computeAllCorHats(tParArray = tParArray,lagCount = lagCount,sample = sample)
      me <- computeMEbyCovHat(tParArray = tParArray,
                              lag = lag,lagCount = lagCount,
                              sample = sample,
                              bandwidth = 0,
                              nonCoverageProbability = nonCoverageProbability,
                              allCorHats = allCor)
      meCenterArray[indexRow,indexCol] <- me[tParArrayCenter]
    }
  }
  fileName <- paste("MEforSampleSize","Alpha",nonCoverageProbability,sep = "_")

  subTitle <- paste(
                    "tPar= ", tParArray[tParArrayCenter],
                    ", Alpha = ", nonCoverageProbability,
                    ", lag = " , lag,
                    sep = "")

  path <- doPath()
  saveCVS(fileName = fileName,path = path,dataToSave = meCenterArray)
  saveJpg(fileName = fileName,path = path)
  matplot(sampleSize,meCenterArray,type = "l",col = 1:nMe)
  title(main = "ME vs SampleSize, fixed alpha and T",sub = subTitle)
  graphics.off()
}
test_that("testing saveMeForSampleSizeFunction",{


  saveMeForSampleSizeFunction()
  cat(" End of test saveMeForSampleSizeFunction", "\n")
  cat("=========================================")
})
