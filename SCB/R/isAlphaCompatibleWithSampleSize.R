#'
#'@export

isAlphaCompatibleWithSampleSize <- function (sampleSize, alpha)
{
  minSampleSize = computeMinSampleSizeForAlpha (alpha)
  isCompatible = sampleSize >= minSampleSize
  return (isCompatible)
}
