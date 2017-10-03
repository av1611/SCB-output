#' @name computeAllBetaLRVHats
#'
#' @aliases computeGammaHat
#'
#' @description Computing a \code{Beta LRV Hat} given certain \code{lag}.
#'
#' @param tParArray a T parameter array, a.k.a. a sleeper parameter, a sequence of numbers {0, 1}.
#' @param lagCount lag count.
#' @param sample a sample.
#' @param kernel a kernel function.
#' @param bandwidth by default, bandwidth = 1.
#'
#' @return A scalar value of a \code{Beta LRV Hat} given certain \code{lag}.
#'
#' @examples
#'

computeAllBetaLRVHats  <- function(tParArray,
                               lagCount,
                               sample,
                               kernel,
                               bandwidth,
                               allCorHats, ...) {
  tParCount <-  length(tParArray)

  allBetaLRVHats <- array(0, dim = c(tParCount, lagCount))
  # array of size tParCount

 for (lagIndex in 1 : lagCount)  {
      allBetaLRVHats[ , lagIndex] <- computeBetaLRVHat(tParArray = tParArray,
                                                      lag = lagIndex,
                                                      sample = sample,
                                                      kernel = kernel,
                                                      bandwidth = bandwidth,
                                                      allCorHats = allCorHats, ...)
    }

  allBetaLRVHats
}
