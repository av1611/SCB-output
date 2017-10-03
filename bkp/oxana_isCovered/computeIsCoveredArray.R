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

computeIsCoveredArray <- function(bandsBrick,
                                  corArray)
{
  replicationCount = dim(bandsBrick)[3]
  mySampleSize = dim(bandsBrick)[2]
  tParCount = dim(bandsBrick)[1]

  isCoveredArray <- array(0, dim = mySampleSize)

  isCoveredBrick <- array(0, dim = c(tParCount, mySampleSize, replicationCount))

  for (replicationIndex in 1:replicationCount)
  {
    band <- bandsBrick[, , replicationIndex]
    isCoveredArray[replicationIndex] <- computeIsCovered(band, corArray)

    isCoveredBrick[, , replicationIndex] <- isCoveredArray[replicationIndex]
  }

  return(isCoveredBrick)
}
