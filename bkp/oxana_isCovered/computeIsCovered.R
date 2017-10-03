#' @title computeIsCovered
#'
#' @name computeIsCovered
#'
#' @aliases
#'
#' @description
#'
#' @param band band
#' @param corArray corArray
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
#' isCovered <- computeIsCovered(band = band, corArray = corArray)


computeIsCovered <- function(band,
                             corArray)
{
  # Array rho_h(i),
  # counter i corresponds to t_i = 0, 1/n, …., 1 - 1/n,
  # horizontal ACF of a sleepers process

  tParCount <- length (corArray)
  isCoveredArrayByT <- array (0, dim = tParCount)

  for (tParIndex in 1 : tParCount)
  {
    upper <- band[, 2]
    lower <- band[, 1]
    if(corArray[tParIndex] <= upper[tParIndex])
    {
      if(corArray[tParIndex] >= lower[tParIndex])
      {
        isCoveredArrayByT[tParIndex] <- 1
      }
      else
      {
        isCoveredArrayByT[tParIndex] <- 0
      }
    }
    else
    {
      isCoveredArrayByT[tParIndex] <- 0
    }
  }

  return(isCoveredArrayByT)
}

