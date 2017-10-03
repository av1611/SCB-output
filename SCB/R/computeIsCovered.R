#' @title computeIsCovered
#'
#' @name computeIsCovered
#'
#' @aliases
#'
#' @description For a given band it checks whether for all t or sleeper values true rho is between upper and lower values. If it is true for all t, it returns 1 or true.
#' If it is not true for at least one t, it returns false or 0.
#'
#' @param band band
#' @param corArray corArray
#' @param fileName fileName
#'
#' @return Vector of {0, 1 values}.
#'
#' @examples
#' tParCount = 10
#' tParArray <- createTParArray(tParCount = tParCount)
#' noise <- createNoise(sampleSize = tParCount, mean = 0, sd = 1)
#' tvMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#' mockTVMA1 <- createTVMA1(sampleSize = tParCount,tvMA1CoefArray = mockTVMA1CoeffArray,noise = mockNoise)
#' corArray <- computeCor(lag = lag, tvMa1CoefArray = tvMA1CoefArray)
#' lag  <- 2
#' band <- createBand(X = rnorm(n = 10, mean = 0, sd = 1), alpha = 0.05)
#' isCovered <- computeIsCovered(band = band, corArray = corArray,fileName=fileName)
#'
#'@export


computeIsCovered <- function(band,
                             corArray,
                             sampleSize,
                             replicationCount,
                             lag,
                             superReplicationCount,
                             nonCoverageProbability,
                             fileName = "")
{
  # Array rho_h(i),
  # counter i corresponds to t_i = 0, 1/n, …., 1 - 1/n,
  # horizontal ACF of a sleepers process
  # saveBand(band = band, corArray = corArray, fileName)

  bandwidth = computeB(sampleSize = sampleSize)
  tParCount <- length (corArray)
  isCoveredArrayByT <- array (0, dim = tParCount)
  saveBand(corArray = corArray,
           band = band,
           sampleSize = sampleSize,
           replicationCount = replicationCount,
           lag = lag,
           superReplicationCount = superReplicationCount,
           nonCoverageProbability = nonCoverageProbability,
           fileName = fileName)

  upperBound <- band[, 2]
  lowerBound <- band[, 1]

  for (tParIndex in 1:tParCount)
  {
    if (corArray[tParIndex] <= upperBound[tParIndex] &
        corArray[tParIndex] >= lowerBound[tParIndex])
      {
        isCoveredArrayByT[tParIndex] <- 1
      }
  }

  isCoveredArrayByT <- prod(isCoveredArrayByT)
}
