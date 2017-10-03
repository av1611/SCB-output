#' @title createDoubleAlphaArray

#' @name createDoubleAlphaArray
#'
#' @aliases createDoubleAlphaArray
#'
#' @description Create double alpha hat array for every alpha in input array
#'
#' @param superReplicationCount how many times we produced alpha hat.
#' @param replicationCount how many bands we need to generate to form a single alphaHat
#' @param sampleSize size of the sample
#' @param lag lag value.
#' @param lagCount how many lags we consider.
#' @param alphaArray a alpha parametr array
#' @param tParArray a T parameter array
#' @param bandwidth a bandwidth value.
#' @return A two dimensional array representing the alpha and many alpha hat for it.
#'
#' @examples
#' myAlphaCount = 4
#'myTParCount  <- 2
#'myTParArray  <- createTParArray(tParCount = myTParCount)
#'mySuperReplicationCount <- 3
#'myReplicationCount <- 6
#'mySampleSize <- 10
#'myLag <- 1
#'myLagCount <- computeLagCount(lag = myLag,sampleSize = mySampleSize)
#'myAlphaArray <- seq(from = 0.05,
#'                  to   = 1 - 1 / myAlphaCount,
#'                   by   = 1 / myAlphaCount)

#'nonCoverageProbabilities <- c(0.2,0.4,0.6,0.8)
#'myBandwidth <- 0.5
#'doubleAlphaArray = createDoubleAlphaArray(
#' superReplicationCount = mySuperReplicationCount,
#'  replicationCount = myReplicationCount,
#' sampleSize = mySampleSize,
#'  alphaArray =myAlphaArray,
#' lag = myLag,
#' lagCount = myLagCount,
#' tParArray = myTParArray,
#'bandwidth = myBandwidth)
#'
#'@export
#'


createDoubleAlphaArray = function(superReplicationCount,
                                  replicationCount,
                                  sampleSize,
                                  lag,
                                  lagCount,
                                  alphaArray,
                                  tParArray,
                                  fileName = "")

{
  isCompatible = isLagCompatibleWithSampleSize(sampleSize = sampleSize, lag = lag)


  alphaCount=length(alphaArray)
  doubleAlphaArray = array(0, dim = c(superReplicationCount, alphaCount))


  for (alphaIndex in 1: alphaCount)
{
    cat("\n Now alpha index =",alphaIndex)
    errorIfNotInputCompatible(alpha = alphaArray[alphaIndex],sampleSize = sampleSize,lag = lag)
    doubleAlphaArray[, alphaIndex] = computeNonCoverageFreqArray(
      superReplicationCount = superReplicationCount,
      replicationCount = replicationCount,
      sampleSize = sampleSize,
      lag = lag,
      lagCount = lagCount,
      tParArray = tParArray,
      nonCoverageProbability = alphaArray[alphaIndex],
      fileName = fileName)
  }


  saveDoubleAplhaHatArray(nonCoverageProbabilities = alphaArray,
                         alphaHats = doubleAlphaArray, sampleSize = sampleSize,
                         lag = lag, replicationCount = replicationCount,
                         superReplicationCount = superReplicationCount,
                         fileName = fileName)

  return(doubleAlphaArray)

}
