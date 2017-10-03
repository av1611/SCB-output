#' @title createSample
#'
#' @name createSample
#'
#' @aliases createX
#'
#' @description Creating \code{sample} given \code{tvMA1CoefArray}tvMA1CoefArray.
#'
#' @param tvMA1CoefArray the psi value.
#'
#' @return The vector of \code{X}
#'
#' @examples
#' createSample(tvMA1CoefArray = seq(from = 0, to = 1, length.out = 10))
#'
#'@export

createSample <- function(sampleSize)
  {
 # tParCount=10 ##Hardcode
  #mockTParArray <- createTParArray(tParCount = tParCount)
  mockTVMA1CoefArray <- createTVMA1CoefArray(sampleSize = sampleSize)
  noise <- createNoise(sampleSize = sampleSize,
                       mean = 0,
                       sd = 1)
  # reserve the space
  sample <- array(0, dim = sampleSize)
  sample[1] <- noise[1]

  for (i in 2:sampleSize)
  {
    sample [i] <- noise[i] + mockTVMA1CoefArray[i] * noise[i - 1]
  }

  sample

}
