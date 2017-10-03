isAlphaCompatibleWithSampleSizeFunction<-function()
{
  cat("\n Testing \'tests_25_isAlphaCompatibleWithSampleSize\' \n")
  alphaArray <- seq(0.1,0.9,by = 0.01)
  requredSampleSize <- numeric(length = length(alphaArray))
  for(alphaIndex in 1:length(alphaArray))
  {
    requredSampleSize[alphaIndex] <- computeMinSampleSizeForAlpha(alpha = alphaArray[alphaIndex])
  }
  df <- cbind(alphaArray,requredSampleSize)
  fileName <- "testIsAlphaCompatibleWithSampleSize"
  saveJpg(fileName = fileName,doPath())
  saveCVS(fileName = fileName,doPath(),df)
  plot(requredSampleSize~alphaArray,type="l")
  graphics.off()
}

test_that("Testing \'isAlphaCompatibleWithSampleSize\'", {
  isAlphaCompatibleWithSampleSizeFunction()
  cat("\n End of test isAlphaCompatibleWithSampleSize", "\n")
  cat("=====================")
})
