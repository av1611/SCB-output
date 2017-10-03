#' @title computeAllCorHats
#'
#' @name computeAllCorHats
#'
#' @aliases computeAllRhoHats
#'
#' @description Computing correlations hat, that is, all \eqn{\rho} hats.
#'
#' @param tParArray a T parameter array, a.k.a. a sleeper parameter, a sequence of numbers {0, 1}.
#' @param lagCount lag.
#' @param sample sample.
#' @param kernel kernel.
#' @param bandwidth bandwidth.
#'
#' @return
#'
#' @examples
#' tParCount = 10
#' tParArray <- createTParArray(tParCount = tParCount)
#' tvMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#' noise <- createNoise(sampleSize = tParCount, mean = 0, sd = 1)
#' mockTVMA1Array <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#' sample <- createSample(model = createTVMA1, tvMA1CoefArray = tvMA1Array, noise = noise)
#' lagCount = tParCount - 1
#' kernel = normalDifferenceKernel
#' allCorHats = computeAllCorHats(tParArray = tParArray, lagCount = lagCount, sample = sample, kernel = kernel,  bandwidth = 1)
#'@export

computeAllCorHats <- function(tParArray,
                              lagCount,
                              sample
                              )

{
  # Double array of sizes sleeperCount and lagCount
  # Inner counter -  tPar
  # Outer counter - lag
  mySampleSize=length(sample)
  bandwidth = computeB(sampleSize = mySampleSize)
  tParCount <- length(tParArray)
  allCorHats <- array(0, dim = c(tParCount, lagCount + 1))
  # colNamesForAllCorHats <- as.character(c(0, seq_len(lagCount)))
  # colnames(allCorHats) <- colNamesForAllCorHats
  # allCorHats <- as.data.frame(allCorHats)
  # allCorHats$`0` <- 1 # replicating 1
  # each tPar element is being called by its index,
  # each lag element is being called as it is
  # (index of each lag is the lag itself)
  for (lagIndex in 0:lagCount) {
    for (tParIndex in seq_len(tParCount)) {
      tParPoint = tParArray[tParIndex]

      # lagIndexToRefer <- eval(parse(text = paste0("allCorHats$`", lagIndex + 1, "`")))
      allCorHats[tParIndex, lagIndex + 1] <- computeCorHat(tParArray = tParPoint,
                                                           lag = lagIndex,
                                                           sample = sample
                                                          )
    }
  }

  allCorHats
  # computeAllCorHats returns one-dimensional array, whereas it should return two-dimensional
}
