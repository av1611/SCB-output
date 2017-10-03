#'
#'@export

isLagCompatibleWithSampleSize <- function (sampleSize, lag)
{
  maxLag = computeMaxLagForSampleSize(sampleSize = sampleSize)
  isCompatible = lag <= maxLag
  return (isCompatible)
}
