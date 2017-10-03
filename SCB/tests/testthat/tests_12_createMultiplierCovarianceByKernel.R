createMultiplierCovarianceByKernelFunction <- function () {
  cat("\n Testing \'createMultiplierCovarianceByKernel\' \n")
  myKernel = normalDifferenceKernel

  mySampleSize = 10
  myBandwidth <- computeB(sampleSize = mySampleSize)
  mockMultiplierCovarianceByKernel <- createMultiplierCovarianceByKernel(kernel = myKernel,
                                                                         bandwidth = myBandwidth,
                                                                         sampleSize = mySampleSize)

  cat ("mockMultiplierCovarianceByKernel[1:3, 1:3] =", mockMultiplierCovarianceByKernel[1:3, 1:3], "\n")
  cat ("dim(mockMultiplierCovarianceByKernel) =", dim(mockMultiplierCovarianceByKernel), "\n")
  expect_that(mockMultiplierCovarianceByKernel, is_a("matrix"))
}

test_that("Testing \'createMultiplierCovarianceByKernel\'", {
  createMultiplierCovarianceByKernelFunction()
  cat(" End of test createMultiplierCovarianceByKernel ","\n")
  cat("=====================")
}
)
