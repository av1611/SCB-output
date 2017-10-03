#'
#'@export

computeMaxLagForSampleSize <- function (sampleSize)
{

  maxLag = sampleSize - computeTermCount(sampleSize) - 1
  return (maxLag)
}
