computeMEForTest <- function(tParArray,
                              lag,
                              lagCount,
                              sample,
                              bandwidth,
                              nonCoverageProbability,
                              allCorHats,
                              C_K = -1.978325,
                              PHI_K_NORMAL_DIFF = 0.4065)
{

  mySampleSize=length(sample)
  betaLRVHat = computeBetaLRVHatForTest(tParArray = tParArray,
                                 lag = lag,
                                 sample = sample,
                                 allCorHats = allCorHats,
                                 bandwidth = bandwidth)


  logSqrt <-  sqrt(-2 * log (bandwidth))
  cFactor <- logSqrt + (C_K - log (log (1 / sqrt (1 - nonCoverageProbability)))) / logSqrt

  sampleSize=length(sample)
  meByCovHat <- cFactor *
    betaLRVHat *
    sqrt(PHI_K_NORMAL_DIFF / (sampleSize * bandwidth))
}
computeBetaLRVHatForTest  <- function(tParArray,
                               lag,
                               sample,
                               bandwidth,
                               allCorHats) {

  mySampleSize=length(sample)
  tParCount = length(tParArray)
  sampleSize = length(sample)
  termCount = floor(2 * sampleSize ^ (4/15)) # aka L

  term=0
  betaLRVHat = array(0, dim = tParCount)
  for (tParIndex in 1 : tParCount)
  {
    for (termIndex in 1 : (termCount))
    {

      term = (2 * allCorHats[tParIndex, lag+1] * allCorHats[tParIndex, termIndex+1] -
                allCorHats[tParIndex, abs(lag - termIndex)+1] -
                allCorHats[tParIndex, lag + termIndex+1]) ^ 2
      betaLRVHat[tParIndex] = betaLRVHat[tParIndex] + term

    }
  }

  betaLRVHat
}
computeAllCorHatsForTest <- function(tParArray,
                              lagCount,
                              sample,
                              bandwidth
)

{

  mySampleSize=length(sample)

  tParCount <- length(tParArray)
  allCorHats <- array(0, dim = c(tParCount, lagCount + 1))
  for (lagIndex in 0:lagCount) {
    for (tParIndex in seq_len(tParCount)) {
      tParPoint = tParArray[tParIndex]
      allCorHats[tParIndex, lagIndex + 1] <- computeCorHatForTest(tParArray = tParPoint,
                                                           lag = lagIndex,
                                                           sample = sample,
                                                           bandwidth = bandwidth)
    }
  }

  allCorHats
  # computeAllCorHats returns one-dimensional array, whereas it should return two-dimensional
}
computeCorHatForTest <- function(tParArray,
                          lag,
                          sample,
                          bandwidth)
{
  kernel=customKernel
  mySampleSize=length(sample)
  myCovariance = computeCovHatForTest(tParArray,
                               lag = lag,
                               sample,
                               bandwidth = bandwidth)


  myVariance = computeCovHatForTest(tParArray,
                             lag = 0,
                             sample,
                             bandwidth = bandwidth)

  corHat = myCovariance / myVariance
}
computeCovHatForTest <- function(tParArray,
                          lag,
                          sample,
                          bandwidth) {
  kernel <- normalDifferenceKernel
  partialSum <- 0
  sampleSize <- length(sample)
  termCountSequence <- seq_len(sampleSize - lag)

  for (termIndex in termCountSequence)
  {
    term <- sample[termIndex] *
      sample[termIndex + lag] *
      kernel((termIndex/sampleSize - tParArray) / bandwidth)

    partialSum <- partialSum + term
  }
  covHat <- partialSum / (sampleSize * bandwidth)
  return (covHat)
}
