computeAllBetaLRVHatsFunction = function () {
  cat("\n Testing \'computeAllBetaLRVHats\' \n")

  myTParCount = 10
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockNoise <- createNoise(sampleSize = myTParCount,
                           mean = 0,
                           sd = 1)
  mockTVMA1Array <- createTVMA1CoefArray(coefFunction = sin,
                                         tParArray = mockTParArray)
  mockSample <- createSample(model = createTVMA1,
                             tvMA1CoefArray = mockTVMA1Array,
                             noise = mockNoise)

  myLagCount = myTParCount - 1
  myKernel = normalDifferenceKernel
  myBandwidth = 1
 # form all rho hats
  # first try fakes
  mockAllCorHats <- computeAllCorHats(tParArray = mockTParArray,
                                    lagCount = myLagCount,
                                    sample = mockSample,
                                    kernel = myKernel,
                                    bandwidth = myBandwidth)


  mockAllBetaLRVHats <- computeAllBetaLRVHats(tParArray = mockTParArray,
                                      lagCount = myLagCount,
                                      sample = mockSample,
                                      kernel = myKernel,
                                      bandwidth = myBandwidth,
                                      allCorHats = mockAllCorHats)

  cat("mockAllBetaLRVHats[1:3] =", mockBetaLRVHat[1:3], "\n")
}

test_that("Testing \'computeAllBetaLRVHats\'", {
  computeAllBetaLRVHatsFunction()
  }
)



