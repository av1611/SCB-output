saveMeForAlphaFunction <- function()
{
  cat ("\n Testing \'tests_74_saveMeForAlpha\'\n")
  sampleSize <- 100
  sample <- createSample(sampleSize = sampleSize)

  nMe <- 5

  nonCoverageProbability <- 0.1

  lag <- 2
  lagCount <- computeLagCount(sampleSize = sampleSize, lag = lag)

  tParCount <- 10
  tParArray <- createTParArray(tParCount = tParCount)

  errorIfNotInputCompatible(alpha = nonCoverageProbability,
                            lag = lag,sampleSize = sampleSize)




  sampleArray <- matrix(nrow = sampleSize,ncol = nMe)
  meArray <- matrix(nrow = length(tParArray),ncol = nMe)


  for(indexCol in 1:nMe)
  {

      sampleArray[,indexCol] <- createSample(sampleSize)

  }
  for(indexCol in 1:nMe)
  {

    allCorHats <- computeAllCorHats(tParArray = tParArray,
                                  lagCount = lagCount,
                                  sample = sampleArray[,indexCol])

    meArray[,indexCol] <- computeMEbyCovHat(tParArray = tParArray,
                            lag = lag,
                            lagCount = myLagCount,
                            sample = sampleArray[,indexCol],
                            nonCoverageProbability = nonCoverageProbability,
                            allCorHats = allCorHats,
                            C_K = -1.978325,
                            # int_sq_der = 0.306951,
                            PHI_K_NORMAL_DIFF = 0.4065)
  }
  fileName <- paste("meArray","Alpha",nonCoverageProbability,sep = "_")

  subTitle <- paste("sampleSize = ", sampleSize,
                    ", tParCount= ", tParCount,
                    ", Alpha = ", nonCoverageProbability,
                    ", lag = " , lag,
                    sep = "")

  path <- doPath()
  saveCVS(fileName = fileName,path = path,dataToSave = meArray)
  saveJpg(fileName = fileName,path = path)
  matplot(tParArray,meArray,type = "l",col = 1:nMe,ylab = " ME",xlab = "Tpar")
  title(main = "ME vs t par, fixed alpha",sub = subTitle)
  dev.off()
}


test_that("testing saveMeForAlpha",{


  saveMeForAlphaFunction()
  cat(" End of test saveMeForAlphaFunction", "\n")
  cat("=========================================")
})
