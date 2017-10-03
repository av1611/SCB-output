computeIsCoveredFunction = function()
{
  cat ("\n Testing \'computeIsCovered\'\n")

  tParCount = 10
  # lagCount = 10
  myLag = 1
  mySampleSize = 10

  myCoefFunction = sin
  mockTParArray <- createTParArray(tParCount = tParCount)

  # mockBand <- createBand(tParArray = mockTParArray,
  #                        lag = myLag,
  #                        lagCount = lagCount,
  #                        bandwidth = 1,
  #                        kernel = normalDifferenceKernel,
  #                        nonCoverageProbability = 0.05)
  lowerBound <- seq(from = 0, to = 1, by = 0.1)
  upperBound <- seq(from = 1, to = 2, by = 0.1)

  mockBand <- cbind(lowerBound, upperBound)

  mockTParArray <- createTParArray(tParCount = tParCount)

  mockTVMA1CoefArray <- createTVMA1CoefArray(coefFunction = myCoefFunction,
                                              sampleSize = mySampleSize)

  mockCor <- computeCor(lag = myLag,
                        tvMa1CoefArray = mockTVMA1CoefArray)
cat("\n Cor",mockCor)
  isCovered <- computeIsCovered(mockBand,
                                mockCor)

  cat("IsCovered:", isCovered, "\n")
}

test_that("Testing computeIsCovered", {computeIsCoveredFunction()})
