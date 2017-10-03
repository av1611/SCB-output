if(!require("lattice")) {
  install.packages("lattice")
  library(lattice)
}


saveMe3DForBandwidthFunction <- function()
{

  workingDirectory<-getwd()
  splitDirectory <- data.frame(strsplit(workingDirectory,"/"))
  tailDirectory <- tail(splitDirectory,1)

  myPath <- "../../../util/meWithBreakedBandwidth.R"

  if(tailDirectory=="SCB")
  {
    myPath <-"../util/meWithBreakedBandwidth.R"
  }
  source(myPath)
  cat ("\n Testing \'tests_76_saveMe3DForBandwidth\'\n")
  sampleSize <- 200
  sample <- createSample(sampleSize = sampleSize)

  tParCount <- 10
  tParArray <- createTParArray(tParCount = tParCount)
  tParArrayCenter <- round(length(tParArray)/2,0)

  alphas <- seq(0.1,0.9,by = 0.1)

  lag <- 1
  lagCount <- computeLagCount(sampleSize = sampleSize, lag = lag)


  bandwidth <- seq(0.1,0.9,by = 0.1)

  meArray <- matrix(nrow = length(bandwidth),ncol = length(alphas))
  for(indexRow in 1:length(bandwidth))
  {
    for(indexCol in 1:length(alphas))
    {
      errorIfNotInputCompatible(alpha = alphas[indexCol],
                                lag = lag,sampleSize = sampleSize)
      allCorHats <- computeAllCorHatsForTest(tParArray = tParArray,
                                      lagCount = lagCount,
                                      sample = sample,
                                      bandwidth = bandwidth[indexRow]
                                      )

      me <- computeMEForTest(tParArray = tParArray,
                                              lag = lag,
                                              lagCount = myLagCount,
                                              sample = sample,
                                              nonCoverageProbability = alphas[indexCol],
                                              allCorHats = allCorHats,
                                              bandwidth = bandwidth[indexRow],
                                              C_K = -1.978325,
                                              # int_sq_der = 0.306951,
                                              PHI_K_NORMAL_DIFF = 0.4065)
      meArray[indexRow,indexCol] <- me[tParArrayCenter]

    }
  }
  fileName <- paste("meArray3DOfBandwidth",sep = "_")

  path <- doPath()
  saveCVS(fileName = fileName,path = path,dataToSave = meArray)
  saveJpg(fileName = fileName,path = path)
  wireframe(meArray,row.values=bandwidth,column.values=alphas,
            xlab="bandwidth",ylab = "NonCoverageProbability",
            scales=list(arrows=FALSE),
            screen = list(z=250,x=-75))
  graphics.off()
}
test_that("testing saveMe3DForBandwidth",
          {
            saveMe3DForBandwidthFunction()
            cat(" End of test saveMe3DForBandwidth", "\n")
            cat("=========================================")
          })
