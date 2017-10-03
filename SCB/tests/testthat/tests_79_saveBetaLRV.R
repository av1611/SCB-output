saveBetaLRVFunction <- function()
{
  cat ("\n Testing \'tests_79_saveBetaLRV\'\n")
  tParCount <- 10
  tParArray <- createTParArray(tParCount)
  sampleSize <- 500
  psi <- customCoefFunction(tParArray)
  betaLrv <-1 + psi^2 / (1 + psi^2)^2
  sample <- createSample(sampleSize)
  lag <- 2
  lagCount <- computeLagCount(sampleSize = sampleSize, lag = lag)
  allCorHat <- computeAllCorHats(tParArray,lagCount,sample)
  betaLrvHat<-computeBetaLRVHat(tParArray = tParArray,lag = lag,
                                 sample = sample,allCorHats =  allCorHat)
  fileName <- "bettaLRVAndBettaLRVHat"
  subTitle <- paste("sampleSize = ", sampleSize,
                    ", tParCount= ", tParCount,
                    ", lag =", lag,
                    sep = "")
  df <- data.frame(betaLrv,betaLrvHat)
  max <- max(df)
  min <- min(df)
  path <- doPath()
  saveCVS(fileName = fileName,path = path,dataToSave = df)
  saveJpg(fileName = fileName,path = path)
  plot(betaLrv ~ tParArray, type = "l",col = "red",ylim = c(min,max))
  lineArray <- c("BetaLRV","BetaLRVHat")
  legend("bottomright",
         title = "LRV",
         lineArray,
         fill = c("red","blue"))
  lines(betaLrvHat~tParArray,col="blue")
  title("True BetaLRV And BetaLRVHat",sub = subTitle)
  graphics.off()
}
test_that("testing saveBetaLRVFunction",
          {
            saveBetaLRVFunction()
            cat(" End of test saveBetaLRVFunction", "\n")
            cat("=========================================")
          })
