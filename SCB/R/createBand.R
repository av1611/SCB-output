#' @title createBand

#' @name createBand
#'
#' @aliases createConfidenceBand
#'
#' @description Create a confidence intervals band.
#'
#' @param lag lag value.
#' @param sample a vector of \code{sample}.
#' @param corArray a vector of \code{corArray}.
#' @param bandwidth a kernel bandwidth value.
#' @param nonCoverageProbability the probability of non-covering the real values with the estimated ones. \code{alpha = 0.05} by default.
#'
#' @return A two dimensional vector representing the lower and the upper bound of CI.
#'
#' @examples
#'  myLag = 1
#'myLagCount = 4
#'myKernel = normalDifferenceKernel
#'myBandwidth = 0.5
#'myNonCoverageProbability = 0.05
#'Start=Sys.time()
#'mySample <- createSample(sampleSize = mySampleSize)
#'  band <- createBand(
#'tParArray = mockTParArray,
#'lag = myLag,
#'lagCount = myLagCount,
#'bandwidth = myBandwidth,
#'kernel = myKernel,
#'sampleSize = mySampleSize,
#'nonCoverageProbability = myNonCoverageProbability)
#'
#'@export

createBand <- function(tParArray,
                       lag,
                       lagCount,
                       sampleSize,
                       nonCoverageProbability)
{
  # get the real correlation computeCor
  # call ME (sample, bandwidth, nonCoverageProbability)
  # upper = corArray + ME
  # upper = corArray + ME

  mySample=createSample(sampleSize = sampleSize)
  mySampleSize=length(mySample)
  bandwidth = computeB(sampleSize = mySampleSize)
  mockallCorHat = computeAllCorHats(
    tParArray = tParArray,
    lagCount = lagCount,
    sample = mySample


  )

  meByCovHat <- computeMEbyCovHat(
    tParArray = tParArray,
    lag = lag,
    lagCount = lagCount,
    sample = mySample,
    nonCoverageProbability = nonCoverageProbability,
    allCorHats = mockallCorHat,
    C_K = -1.978325,
    PHI_K_NORMAL_DIFF = 0.4065)

  corHat = computeCorHat(
    tParArray = tParArray,
    lag = lag,
    sample = mySample

  )


  lowerBound <- corHat - meByCovHat
  upperBound <- corHat + meByCovHat
  band <- cbind(lowerBound, upperBound)
  # saveBand(corArray = corArray,band =band,fileName = myfileName )
  return(band)



}
