if(!require("lattice")) {
  install.packages("lattice")
  library(lattice)
}


computeLagCountFunction <- function()
{

  cat("\n Testing \'tests_22_computeLagCount\' \n")
  maxGraph=1000
  step=100
  sampleSizes <- seq(1,1000,by=100)
  lags <- seq(1,10,by = 1)
  alpha <- 0.1

  lagCounts <- matrix(nrow = length(sampleSizes),ncol = length(lags))

  for(i in 1:length(sampleSizes))
  {
    for(j in 1:length(sampleSizes))
    {
      lagCounts[i,j]<- computeLagCount(sampleSize = sampleSizes[i],lag = lags[j])
    }
  }

  fileName <- "LagCount"

  df <- data.frame(lagCounts)

  saveCVS(fileName = fileName, doPath(), df)

  saveJpg(fileName = fileName, doPath())



  wireframe(lagCounts,row.values = sampleSizes,column.values = lags,
            xlab="sampleSizes",ylab="Lags",
            screen = list(z=-300,x=-75),
            scales=list(arrows=FALSE))

  graphics.off()

}

test_that("Testing computeLagCount", {
  computeLagCountFunction()
  cat("\n End of test computeLagCount", "\n")
  cat("=====================")
})

