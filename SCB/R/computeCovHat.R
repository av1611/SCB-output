#' @title computeCovHat
#'
#' @name computeCovHat
#'
#' @aliases
#'
#' @description
#'
#' @param tParArray a T parameter, a.k.a. a sleeper parameter, a sequence of numbers {0, 1}.
#' @param lag a lag value.
#' @param sample a sample.
#' @param kernel a kernel function.
#' @param bandwidth by default, bandwidth = 1.
#'
#' @return
#'
#' @examples
#' tParArray <- createTParArray(tParCount = tParCount)
#' sample <- createSample(model = createMa1, tvMA1CoefArray = tvMA1CoefArray, noise = noise)
#' kernel  <- normalDiferencekernel
#' bandwith <- 1
#' lag <- 2
#'
#'@export

computeCovHat <- function(tParArray,
                          lag,
                          sample) {
  sampleSize <- length(sample)
  bandwidth <- computeB(sampleSize = sampleSize)
  partialSum <- 0
  termCountSequence <- seq_len(sampleSize - lag)

    for (termIndex in termCountSequence)
    {
     term <- sample[termIndex] *
                sample[termIndex + lag] *
                customKernel((termIndex/sampleSize - tParArray) / bandwidth)

      partialSum <- partialSum + term
    }




  covHat <- partialSum / (sampleSize * bandwidth)
  return (covHat)
}



