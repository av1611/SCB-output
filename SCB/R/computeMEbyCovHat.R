#' @title computeMEbyCovHat
#'
#' @name computeMEbyCovHat
#'
#' @aliases
#'
#' @descriptionHeart Heart of the whole thesis.Call GammaHatCompute CK.Compute psiK
#'
#' @param
#' tParArray,
#' lag,
#' lagCount,
#' sample,
#' bandwidth,
#' nonCoverageProbability,
#' allCorHats,
#' C_K = -1.978325,
#' PHI_K_NORMAL_DIFF = 0.4065
#' @return Margin of error
#'
#' @examples
#' tParCount = 10
#' tParArray <- createTParArray(tParCount = tParCount)
#' noise <- createNoise(sampleSize = tParCount, mean = 0, sd = 1)
#' tvMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#' sample <- createSample(model = createMa1, tvMA1CoefArray = tvMA1CoefArray, noise = noise)
#' bandwith <- 0.5
#' lag <- 2
#' nonCoverageProbability <- 0.05
#' meByCovHat <- computeMEbyCovHats(lag = lag, sample = sample, bandwith = bandwidth, nonCoverageProbability = nonCoverageProbability
#'
#'@export

computeMEbyCovHat <- function(tParArray,
                              lag,
                              lagCount,
                              sample,
                              bandwidth,
                              nonCoverageProbability,
                              allCorHats,
                              C_K = -1.978325,
                              PHI_K_NORMAL_DIFF = 0.4065)
  {
  # Array rho_h(i),
  # counter i corresponds to t_i = 0, 1/n, …., 1 - 1/n,
  # horizontal ACF of a sleepers process
  #lengthSleeper <- length(sleeper)
  # do we need this?
  # rhoHI <- rep(0, 10)
  # mylag=2
  # int_sq_der = 0.306951 is now hidden in C_K

  mySampleSize=length(sample)
  bandwidth = computeB(sampleSize = mySampleSize)
  betaLRVHat = computeBetaLRVHat(tParArray = tParArray,
                                 lag = lag,
                                 sample = sample,
                                 allCorHats = allCorHats)


  logSqrt <-  sqrt(-2 * log (bandwidth))
  cFactor <- logSqrt + (C_K - log (log (1 / sqrt (1 - nonCoverageProbability)))) / logSqrt

  sampleSize=length(sample)
  meByCovHat <- cFactor *
                betaLRVHat *
                sqrt(PHI_K_NORMAL_DIFF / (sampleSize * bandwidth))
  }
