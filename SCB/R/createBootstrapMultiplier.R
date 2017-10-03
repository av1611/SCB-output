#' @title createBootstrapMultiplier
#'
#' @name createBootstrapMultiplier
#'
#' @aliases createBM
#'
#' @description Create bootstrap W series using Cholesky decomposition
#'
#' @param kernel kernel.
#' @param bandwidth bandwidth = 1, by default
#' @param sampleSize the sample size.
#'
#' @return The covariance matrix.
#'
#' @examples
#' createBootstrapMultiplier(kernel, bandwidth = 1, sampleSize = 10)
#'
#'@export

createBootstrapMultiplier <- function(kernel = customKernel,
                                      bandwidth = 1,
                                      sampleSize) {
  noise = createNoise(sampleSize, 0, 1)
  covar = createMultiplierCovarianceByKernel(kernel,
                                             bandwidth,
                                             sampleSize)
  W = chol(covar) %*% noise
  W
}
