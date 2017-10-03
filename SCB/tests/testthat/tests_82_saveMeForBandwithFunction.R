saveMeForBandwithFunction <- function()
{
  cat ("\n Testing \'tests_82_saveMeForBandwith\'\n")
  workingDirectory<-getwd()
  splitDirectory <- data.frame(strsplit(workingDirectory,"/"))
  tailDirectory <- tail(splitDirectory,1)

  myPath <- "../../../util/meWithBreakedBandwidth.R"

  if(tailDirectory=="SCB")
  {
    myPath <-"../util/meWithBreakedBandwidth.R"
  }
  source(myPath)
  sampleSize <- 500
  bandwith <- seq(0.1,0.9,0.2)
  tParCount <- 10
  tParArray <- createTParArray(tParCount)

  nonCoverageProbability <- 0.2
  lag <- 1
  errorIfNotInputCompatible(alpha = nonCoverageProbability,lag = lag,sampleSize = sampleSize)
  lagCount <- computeLagCount(sampleSize = sampleSize,lag = lag)
  sample <- createSample(sampleSize)
  meArray <- matrix(0,nrow = length(tParArray),ncol = length(bandwith))
  for(index in 1:length(bandwith))
  {
    tempBandwith <- bandwith[index]
    corArray <- computeAllCorHatsForTest(tParArray = tParArray,
                                         lagCount = lagCount,
                                         sample = sample,
                                         bandwidth = tempBandwith)
    meArray[,index] <- computeMEForTest(tParArray = tParArray,
                                lag = lag,
                                lagCount = lagCount,
                                sample = sample,
                                bandwidth = tempBandwith,
                                nonCoverageProbability = nonCoverageProbability,
                                allCorHats = corArray)
  }
  fileName <- paste("MEforBandWith","Alpha",nonCoverageProbability,"SS",sampleSize,sep = "_")
  xlab <- "TparArray"
  subTitle <- paste(xlab,
    "\ntParCount= ", tParCount,
    ", Alpha = ", nonCoverageProbability,
    ", lag = " , lag,
    ", sampleS = ", sampleSize,
    "\nAlpha, sample,SampleSize fixed",
    sep = "")

  path <- doPath()
  saveCVS(fileName = fileName,path = path,dataToSave = meArray)
  saveJpg(fileName = fileName,path = path)
  matplot(tParArray,meArray,type = "l",col = 1:length(bandwith),
          xlab = "")
  title(main = "ME(t) for different bandwidth",sub = subTitle)
  legend("topright",title = "bandwith",legend = bandwith,
         fill = 1:length(bandwith))
  graphics.off()
}
test_that("testing saveMeForBandwithFunction",{


  saveMeForBandwithFunction()
  cat(" End of test saveMeForBandwithFunction", "\n")
  cat("=========================================")
})

