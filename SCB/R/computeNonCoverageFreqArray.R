#' @title computeNonCoverageFreqDist
#'
#' @name computeNonCoverageProbDist
#'
#' @aliases
#'
#' @description in a cycle for each alpha in array call NonCoverageFreqArray
#'
#'
#' @param superReplicationCount superReplicationCount
#' @param replicationCount a replication count.
#' @param sampleSize a sample size.
#' @param lag a lag value, true lag is lag - 1.
#' @param tParArray a T parameter array
#' @param bandwidth a bandwidth value.
#' @param nonCoverageProbability probability of non-coverage.
#'
#' @return array of alpha hat
#'
#' @examples
#' myTParCount  <-2
#' myTParArray  <- createTParArray(tParCount = myTParCount)
#'mySuperReplicationCount <- 4
#'myReplicationCount <- 4
#'mySampleSize <- 50
#'myLag <- 1
#'myLagCount <- 4
#'myNonCoverageProbability <- 0.05
#'nonCoverageFreqArray <- computeNonCoverageFreqArray(
#'superReplicationCount = mySuperReplicationCount,
#' replicationCount = myReplicationCount,
#'  sampleSize = mySampleSize,
#' lag = myLag,
#'lagCount = myLagCount,
#'  tParArray = myTParArray,
#' bandwidth = myBandwidth,
#' nonCoverageProbability = myNonCoverageProbability)
#'
#'@export



computeNonCoverageFreqArray  <- function(superReplicationCount,
                                        replicationCount,
                                        sampleSize,
                                        lag,
                                        lagCount,
                                        tParArray,
                                        nonCoverageProbability,
                                        fileName = "") {

  bandwidth = computeB(sampleSize = sampleSize)
  alphaHatArray <- array(0, dim = superReplicationCount)

  corArray <- computeCor(lag = lag,tParArray = tParArray)

  for (superIndex in 1:superReplicationCount) {

    cat("\n Now super index =",superIndex)

    alphaHatArray[superIndex] <- computeNonCoverageFreq(replicationCount=replicationCount,
                                                        sampleSize = sampleSize,
                                                        lag = lag,
                                                        tParArray = tParArray,
                                                        corArray = corArray,
                                                        lagCount = lagCount,
                                                        nonCoverageProbability = nonCoverageProbability,
                                                        superReplicationCount = superReplicationCount,
                                                        fileName = fileName)

  }






   saveNonCoverageFreqArray(nonCoverageProbability = nonCoverageProbability,
                            alphaHatArray = alphaHatArray,
                            sampleSize = sampleSize, bandwidth = bandwidth,
                            replicationCount = replicationCount, lag = lag,
                            superReplicationCount = superReplicationCount,
                            fileName = fileName)

  alphaHatArray



}

