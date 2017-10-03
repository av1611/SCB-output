#' @title computeCorHat
#'
#' @name computeCorHat
#'
#' @aliases computeRhoHat
#'
#' @description Computing a ratio of \code{correlation hat} given certain \code{lag} over \code{correlation hat} of \code{lag = 0}.
#'
#' @param tParArray a T parameter, a.k.a. a sleeper parameter, a sequence of numbers {0, 1}.
#' @param lag a lag value.
#' @param sample a sample.
#' @param bandwidth by default, bandwidth = 1.
#'
#' @return A scalar value of a ratio of \code{correlation hat} given certain \code{lag} over \code{correlation hat} of \code{lag = 0}
#'
#' @examples
#' tParCount = 10
#' tParArray <- createTParArray(tParCount = tParCount)
#' noise <- createNoise(sampleSize = tParCount, mean = 0, sd = 1)
#' tvMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#' sample <- createSample(model = createMa1, tvMA1CoefArray = tvMA1CoefArray, noise = noise)
#' lag = 2
#' corHat = computeCorHat(tParArray = tParArray, lag = lag, sample = sample, bandwidth = bandwith)
#'
#'@export


computeCorHat <- function(tParArray,
                          lag,
                          sample)
{
  mySampleSize=length(sample)
  bandwidth = computeB(sampleSize = mySampleSize)
  myCovariance = computeCovHat(tParArray,
                               lag = lag,
                               sample)
  # myCovlen=length(myCovariance)
  # cat("\n len myCov= ",myCovlen)
  # cat("\n myCovar= ",myCovariance)
  myVariance = computeCovHat(tParArray,
                             lag = 0,
                             sample)
  # myVarlen=length(myVariance)
  # cat("\n len myVar= ",myVarlen)
  # cat("\n myVar= ",myVariance)
  corHat = myCovariance / myVariance
}
