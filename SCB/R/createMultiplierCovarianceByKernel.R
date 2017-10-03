#' @title createMultiplierCovarianceByKernel
#'
#' @name createMultiplierCovarianceByKernel
#'
#' @aliases createCovarianceByKernel
#'
#' @description Generate covariance function given kernel. Is being called from within \code{createBootstrapMultiplier} function.
#'
#' @param kernel kernel.
#' @param bandwidth bandwidth = 1, by default
#' @param sampleSize the sample size.
#'
#' @return The covariance matrix.
#'
#' @examples
#' createMultiplierCovarianceByKernel(kernel, bandwidth = 1, sampleSize = 10)
#'
#'@export

createMultiplierCovarianceByKernel <- function(kernel = customKernel,
                                               bandwidth = 1,
                                               sampleSize) {
  covar = diag(kernel(0/bandwidth)/2, sampleSize)

  for (i in 1:sampleSize - 1) {
    for (j in (i + 1):sampleSize) {
      covar[i, j] = kernel(abs(j - i)/bandwidth)
      }
  }

  covar = covar + t(covar)
}
