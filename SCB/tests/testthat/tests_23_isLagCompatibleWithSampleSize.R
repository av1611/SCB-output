isLagCompatibleWithSampleSizeFunction <- function(maxGraph = 100,step = 10)
{
  cat("\n Testing \'tests_23_isLagCompatibleWithSampleSize\' \n")

  sampleSizes <- seq(0,maxGraph,by=step)
  lags <- seq(0,maxGraph,by = step)

  compatibles <- matrix(nrow = length(sampleSizes),ncol = length(lags))
  cocoses <- numeric(length = length(sampleSizes))

  for(i in 1:length(lags))
  {
    for(j in 1:length(sampleSizes))
    {
      compatibles[i,j] <- isLagCompatibleWithSampleSize(sampleSize = sampleSizes[j], lag = lags[i])

      if(isLagCompatibleWithSampleSize(sampleSize = sampleSizes[j],lag = lags[i]))
      {
        cocoses[i] <- j
        break
      }
    }
  }

  fileName <- "IsinPutCompatible"

  xMaxLim <- floor(maxGraph * 0.8)
  yMaxLim <- max(cocoses)

  saveJpg(fileName = fileName, doPath())
  saveCVS(fileName = fileName, doPath(), compatibles)

  plot(cocoses ~ lags,type = "l",xlim=c(0,xMaxLim),ylim=c(0,yMaxLim),
       col="red",ylab = "sampleSize",xlab = "lag")
  title("Is input Compatible")
  legend("bottomright",c("compatibles"),fill = c("red"))
  #lines(cocoses,lags,col="red")

  graphics.off()
}
test_that("Testing \'isLagCompatibleWithSampleSize\'", {
  isLagCompatibleWithSampleSizeFunction()
  cat("\n End of test isLagCompatibleWithSampleSize", "\n")
  cat("=====================")
})

