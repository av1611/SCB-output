computeIsCoveredArrayFunction <- function()
{
  cat("\n Testing \'computeIsCoveredArray\'\n")

  myLag = 1
  tParCount = 2
  mySampleSize = 2
  replicationCount = 3
  myCoefFunction = sin

  mockTParArray <- createTParArray(tParCount = tParCount)

  mockTVMA1CoefArray <-
    createTVMA1CoefArray(coefFunction = myCoefFunction,
                         sampleSize = mySampleSize)

  mockCorArray <- computeCor(lag = myLag,
                             tvMa1CoefArray = mockTVMA1CoefArray)

  lowerBound <- seq(from = 0, to = 1, by = 0.1)
  upperBound <- seq(from = 1, to = 2, by = 0.1)
  mockBand <- cbind(lowerBound, upperBound)

  mockBandsBrick <-
    array(mockBand, dim = c(tParCount, mySampleSize, replicationCount))

  # mockBandsBrick <- createBandsBrick(sampleSize = tParCount,
  #                                    replicationCount = tParCount)

  isCoveredArray <-
    computeIsCoveredArray(bandsBrick = mockBandsBrick,
                          corArray = mockCorArray)

  cat("isCoveredArray:", isCoveredArray, "\n")
}

test_that("Testing \'isCoveredArray\'", {
  computeIsCoveredArrayFunction()
})
