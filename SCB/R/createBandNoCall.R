#'
#'@export

createBandFunctionNoCall <- function()
{
  mySampleSize=5
  myTParCount = 10
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockTVMA1Array <- createTVMA1CoefArray(coefFunction = customCoefFunction,
                                         sampleSize = mySampleSize)

  mockSample <- createSample(sampleSize = myTParCount)

  myLag = 3
  myLagCount = myTParCount - 1
  myBandwidth = 1
  # form all rho hats
  # first try fakes
  mockAllCorHats <- computeAllCorHats(tParArray = mockTParArray,
                                      lagCount = myLagCount,
                                      sample = mockSample,
                                      bandwidth = 0.5)
cat("\n allCorHats",str(mockAllCorHats))

meByCovHat <- computeMEbyCovHat(tParArray = mockTParArray,
                                lag = myLag,
                                lagCount = myLagCount,
                                sample = mockSample,
                                bandwidth = 0.5,
                                nonCoverageProbability = 0.05,
                                allCorHats = mockAllCorHats,
                                C_K = -1.978325,
                                PHI_K_NORMAL_DIFF = 0.4065)
cat("\n ME",meByCovHat)

corHat <- computeCorHat(tParArray = mockTParArray,
                        lag = myLag,
                        sample = mockSample,
                        bandwidth = 0.5)

cat("\n corHat",corHat)

lowerBound <- corHat - meByCovHat
upperBound <- corHat + meByCovHat
band <- cbind(lowerBound, upperBound)
cat("\n band ",band)
  # # expect_that(dim(mockBand)[1], equals(2))  # the number of rows
  # # expect_that(dim(mockBand)[2], equals(10)) # the number of cols
}


