#' @title computeNonCoverageFreq
#'
#' @name computeNonCoverageFreq
#'
#' @aliases
#'
#' @description For each band in bandsArray it computes isCovered and ads up all those zeros and ones. Divides this sum by replicationCount = length or dim of bandsArray
#'
#' @param replicationCount a replication count.
#' @param sampleSize a sample size.
#' @param lag a lag value, true lag is lag - 1.
#' @param tParArray a T parameter array
#' @param corArray a true correlation array
#' @param bandwidth a bandwidth value.
#' @param nonCoverageProbability probability of non-coverage.
#'
#' @returns A scalar value of non-coverage frequency, zero count over replication count.
#'
#' @examples
#'myTParCount  <-2
#'mockTParArray  <- createTParArray(tParCount = myTParCount)
#'myReplicationCount <- 20
#'mySampleSize <-50
#'myLag <- 2
#'myLagCount <- 4
#'trueCorArray <- computeCor(lag = myLag,coefFunction = sin,tParArray = mockTParArray)
#'myBandwidth <- 0.5
#'myNonCoverageProbability <- 0.05
#'bandsBrick = createBandsBrick(sampleSize = mySampleSize,
#'                                tParArray  = mockTParArray,
#'                                lag        = myLag,
#'                                lagCount   = myLagCount,
#'                                bandwidth  = myBandwidth,
#'                              nonCoverageProbability = myNonCoverageProbability,
#'                               replicationCount       = myReplicationCount)
#'
#'@export




computeNonCoverageFreq  <- function(replicationCount,
                                    sampleSize,
                                    lag,
                                    tParArray,
                                    corArray,
                                    lagCount,
                                    nonCoverageProbability = 0.05,
                                    superReplicationCount,
                                    fileName = "") {
  bandwidth = computeB(sampleSize = sampleSize)
  errorIfNotInputCompatible(alpha = nonCoverageProbability,sampleSize=sampleSize,lag=lag)
  bandsBrick = createBandsBrick(sampleSize = sampleSize,
                                tParArray  = tParArray,
                                lag        = lag,
                                lagCount   = lagCount,
                                bandwidth  = bandwidth,
                                nonCoverageProbability = nonCoverageProbability,
                                replicationCount       = replicationCount)

  isCoveredArray <- computeIsCoveredArray(bandsBrick,
                                          corArray = corArray,tParCount = length(tParArray),
                                          sampleSize = sampleSize,
                                          replicationCount = replicationCount,
                                          lag = lag,
                                          superReplicationCount = superReplicationCount,
                                          nonCoverageProbability = nonCoverageProbability,
                                          fileName = fileName)

  zeroCount = 0
  replicationCount = dim(bandsBrick)[1]

  for (i in 1:length(isCoveredArray)) {
    if (isCoveredArray[i] == 0) {
      zeroCount <- zeroCount + 1
    }
  }

  nonCoverageFreq <- zeroCount / replicationCount

}
