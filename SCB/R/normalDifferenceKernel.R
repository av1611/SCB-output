#' @title normalDifferenceKernel
#'
#' @name normalDifferenceKernel
#'
#' @aliases kernelGaussian4
#'
#' @description A manual implementation of the 4th Gaussian kernel.
#'
#' A manual implementation of the 4th Gaussian kernel.
#'
#' @param u the length.
#'
#' @return A scalar kernel
#'
#' @examples
#' normalDifferenceKernel(10)
#'
#'@export

normalDifferenceKernel  <-  function(u)  {
  #4th Gaussian
  2 * dnorm(u) - dnorm(u / sqrt(2)) / sqrt(2)
}

