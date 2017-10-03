errorIfNotInputCompatibleFunction=function()
{
  cat("\n Testing \'tests_06_errorIfNotInputCompatible\'")
  mySampleSize=40
  myLag= 0
  myAlpha <- 0.3
  errorIfNotInputCompatible(alpha = myAlpha,sampleSize = mySampleSize,lag=myLag)

}
test_that("Testing \'errorIfNotInputCompatibleFunction\'", {
  errorIfNotInputCompatibleFunction()
  cat("\n End of test errorIfNotInputCompatible", "\n")
  cat("=====================")
})
