saveBandwithFunction <- function()
{
  cat("\n Testing \'tests_24_saveBandwith\' \n")
  sampleSize <- 1:500

  bandwidth <- numeric(length = length(sampleSize))

  for(index in sampleSize)
  {
    bandwidth[index] <- computeB(sampleSize[index])
  }
  path <- doPath()

  saveJpg(fileName = "Bandwith",path = path)
  plot(bandwidth ~ sampleSize, type = "l", main = "Bandwith of samplesize")
  graphics.off()
}

test_that("Testing \'saveBandwithFunction\'", {
  saveBandwithFunction()
  cat("\n End of test saveBandwithFunction", "\n")
  cat("=====================")
})
