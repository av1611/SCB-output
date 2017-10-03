createBootstrapMultiplierFunction <- function () {
  cat("\n Testing \'createBootstrapMultiplier\' \n")

  myKernel = normalDifferenceKernel

  mySampleSize = 10
  myBandwidth <- computeB(sampleSize = mySampleSize)
  mockBootstrapMultiplier <- createBootstrapMultiplier(kernel = myKernel,
                                                       bandwidth = myBandwidth,
                                                       sampleSize = mySampleSize)

  cat ("mockBootstrapMultiplier[1:3] =", mockBootstrapMultiplier[1:3], "\n")
  cat ("dim(mockBootstrapMultiplier) =", dim(mockBootstrapMultiplier), "\n")
  expect_that(mockBootstrapMultiplier, is_a("matrix"))
}

test_that("Testing \'createBootstrapMultiplier\'", {
  createBootstrapMultiplierFunction()
  cat(" End of test createBootstrapMultiplier","\n")
  cat("=====================")
}
)
