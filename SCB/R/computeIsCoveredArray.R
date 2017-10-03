#' @title computeIsCoveredArray
#'
#' @name computeIsCoveredArray
#'
#' @aliases
#'
#' @description in a cycle for each band in a brick call computeIsCovered
#'
#' @param bandsBrick bandsBrick
#' @param corArray corArray
#'
#' @return Brick of vectors with values {0, 1}.
#'
#' @examples
#'mySampleSize=5
#'myTParCount = 5
#'mockTParArray <- createTParArray(tParCount = myTParCount)
# may be different
#'mockTVMA1Array <- createTVMA1CoefArray(coefFunction = sin,sampleSize = mySampleSize)
#'myReplicationCount=2
#'myLag = 1
#'myLagCount = 3
#'myKernel = normalDifferenceKernel
#'myBandwidth = 0.5
#'myNonCoverageProbability=0.05
#'fileName <- paste("ss", mySampleSize, "l", myLag, "bandW", myBandwidth, "alpha"
#'                , myNonCoverageProbability, sep = "_")
#'bandsBrick=createBandsBrick(tParArray=mockTParArray,
#'                            lag=myLag,
#'                           lagCount=myLagCount,
#'                            bandwidth= myBandwidth,
#'                            kernel = normalDifferenceKernel,
#'                            sampleSize=mySampleSize,
#'                            nonCoverageProbability = myNonCoverageProbability,
#'                            replicationCount=myReplicationCount)
#'mockCorArray <- computeCor(lag = myLag,coefFunction = sin,tParArray = mockTParArray)
#'isCoveredArray <-
#' computeIsCoveredArray(bandsBrick = bandsBrick,
#'                        corArray = mockCorArray,
#'                       fileName = fileName)
#'
#'@export

computeIsCoveredArray <- function(bandsBrick,
                                  corArray,
                                  tParCount,
                                  sampleSize, bandwidth,
                                  replicationCount,
                                  lag,
                                  superReplicationCount,
                                  nonCoverageProbability,
                                  fileName = "")
{
  # replicationCount = dim(bandsBrick)[1]
  # mySampleSize = dim(bandsBrick)[2]
  # tParCount = dim(bandsBrick)[3]

  isCoveredArray <- array(0, dim = replicationCount)

  isCoveredBrick <-
    array(0, dim = c(tParCount, sampleSize, replicationCount))

  for (replicationIndex in 1:replicationCount)
  {
    band <- bandsBrick[replicationIndex, , ]

    isCoveredArray[replicationIndex] <-
      computeIsCovered(band,
                       corArray = corArray,
                       sampleSize = sampleSize,
                       replicationCount = replicationCount,
                       lag = lag,
                       superReplicationCount = superReplicationCount,
                       nonCoverageProbability = nonCoverageProbability,
                       fileName = fileName)

  }

  isCoveredArray
}
