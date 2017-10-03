saveCorAndCorHatFunction <- function()
{
  cat ("\n Testing \'tests_80_saveCorAndCorHat\'\n")
  sampleSize <- 1000000
  tParCount <- 10
  tParArray <- createTParArray(tParCount = tParCount)


  lag <- 1

  nCorHat <- 5
  corArray <- computeCor(lag,tParArray = tParArray)
  corHatArray <- matrix(0,nrow=length(tParArray),ncol = nCorHat)
  for(index in 1:nCorHat)
  {
    sample <- createSample(sampleSize = sampleSize)
    bandWidth <- computeB(sampleSize = sampleSize)
    corHatArray[,index] <- computeCorHat(tParArray = tParArray,
                                         lag = lag,
                                         sample = sample)
  }
  subTitle <- paste("sampleSize = ", sampleSize,
                    ", tParCount= ", tParCount,
                    ", lag =", lag,
                    sep = "")
  path <- doPath()
  fileName<-"CorAndCorHat"
  tempFindMax <- cbind(corHatArray,corArray)
  saveCVS(fileName = fileName,path = path,dataToSave = tempFindMax)
  max <- max(tempFindMax)
  min <- min(tempFindMax)
  saveJpg(fileName = fileName,path = path)
  matplot(tParArray,corHatArray,type = "l",col = "blue",ylim = c(min,max))
  lines(tParArray,y=corArray,col="red")
  title("Cor and CorHat",sub = subTitle)
  graphics.off()
}
test_that("testing saveCorAndCorHatFunction",
          {
            saveCorAndCorHatFunction()
            cat(" End of test saveCorAndCorHatFunction", "\n")
            cat("=========================================")
          })

