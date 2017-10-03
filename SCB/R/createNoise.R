#' @title createNoise
#'
#' @name createNoise
#'
#' @aliases createZ
#'
#' @description Creating a normally distributed vector of \code{noise} given \code{sampleSize} (length), \code{mean} and \code{sd}.
#'
#' @param sampleSize time, a.k.a. \code{N} of the sample of \code{X}.
#' @param mean the mean.
#' @param sd standard deviation.
#'
#' @return The vector of numbers.
#'
#' @examples
#' createNoise(sampleSize = 10, mean = 0, sd =  1)
#' createNoise(sampleSize = 20, mean = 5, sd =  10)
#'
#'@export

createNoise <- function(sampleSize,
                        mean,
                        sd) {
  noise <- rnorm(n = sampleSize,
                 mean = mean,
                 sd = sd)
}
