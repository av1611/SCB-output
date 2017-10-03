#'
#'@export

computeLagCount=function (sampleSize, lag)
{
  myTermCount = computeTermCount (sampleSize)
  return (myTermCount + lag)
}
