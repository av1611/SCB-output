#' @title computeBetaLRVHat
#'
#' @name computeBetaLRVHat
#'
#' @aliases computeGammaHat
#'
#' @description Computing a \code{betaLRVHat} given certain \code{lag}.
#'
#' @param tParArray a T parameter array, a.k.a. a sleeper parameter, a sequence of numbers {0, 1}.
#' @param lag a lag value.
#' @param sample a sample.
#' @param bandwidth by default, bandwidth = 1.
#'
#' @return A scalar value of a \code{Beta LRV Hat} given certain \code{lag}.
#'
#' @examples
#' tParCount = 10
#' tParArray <- createTParArray(tParCount = tParCount)
#' noise <- createNoise(sampleSize = tParCount, mean = 0, sd = 1)
#' tvMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#' sample <- createSample(model = createMa1, tvMA1CoefArray = tvMA1CoefArray, noise = noise)
#' bandwith = 1
#' lag = 2
#' allCorHats = computeAllCorHats(tParArray = tParArray, lagCount = tParCount - 1, sample = sample, bandwidth = bandwidth)
#' betaLRVHat = computeBetaLRVHat(tParArray = tParArray, lag = lag, sample = sample, bandwidth = bandwidth, allCorHats = allCorHats)
#'
#'@export

computeBetaLRVHat  <- function(tParArray,
                               lag,
                               sample,
                               allCorHats) {
  mySampleSize=length(sample)
  bandwidth = computeB(sampleSize = mySampleSize)
  tParCount = length(tParArray)
  sampleSize = length(sample)
  #cat("\nSampleSize",sampleSize)
  termCount = floor(2 * sampleSize ^ (4/15)) # aka L
  #cat("\ntermCount",termCount)
  # check sizes
  # if (lag > sampleSize - 1 - termCount) {
  #   cat ("after if with stop \n")
  #   cat ("sampleSize = ",sampleSize, "\n")
  #   cat ("termCount = ",termCount, "\n")
  #   cat ("lag = ",lag, "\n")
  #   cat("wrong relation between sampleSize, termCount, lag! \n")
  #   # stop("wrong relation between sampleSize, termCount, lag!")
  # }

  # if (! lag > termCount) {
  #   cat ("\n after if with stop \n")
  #   stop("wrong relation between sampleSize, termCount, lag!")
  # }
  term=0
  betaLRVHat = array(0, dim = tParCount)
  # array of size tParCount
  #xTest <- termCount - lag -1
  #cat("\n termCount - lag - 1: ",xTest)
  for (tParIndex in 1 : tParCount)
  {
    for (termIndex in 1 : (termCount))
    {

      term = (2 * allCorHats[tParIndex, lag+1] * allCorHats[tParIndex, termIndex+1] -
                allCorHats[tParIndex, abs(lag - termIndex)+1] -
                allCorHats[tParIndex, lag + termIndex+1]) ^ 2
      # minuendFactor1 = 2 * allCorHats[tParIndex, lag+1]
      # minuendFactor2 = allCorHats[tParIndex, termIndex+1]
      # minuend = minuendFactor1 * minuendFactor2
      # subtrahend1 = allCorHats[tParIndex, abs(lag - termIndex)+1]
      # lagIndex = lag + termIndex+1
      # cat ("\n lagIndex = ", lagIndex, ", while max = ", dim (allCorHats)[2])


      # subtrahend2 = allCorHats[tParIndex, lagIndex]
      # term = (minuend - subtrahend1 - subtrahend2)^2
      # cat("\n minuendFactor1 = ", minuendFactor1)
      # cat("\n minuendFactor2 = ", minuendFactor2)
      # cat("\n minuend= ",minuend)
      # cat("\n tParIndex= ",tParIndex)
      # cat("\n termIndex+1= ",termIndex)
      # cat("\n subtrahend1= ",subtrahend1)
      # cat("\n subtrahend2= ",subtrahend2)
      # cat("\n term= ",term)
      # cat("\n termCount= ",termCount)


      betaLRVHat[tParIndex] = betaLRVHat[tParIndex] + term

    }
  }

  betaLRVHat
}
