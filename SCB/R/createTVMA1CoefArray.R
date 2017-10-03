#' @title createTVMA1CoefArray
#'
#' @name createTVMA1CoefArray
#'
#' @aliases createPsi
#'
#' @description Creating \code{ma1CoefArray}, that is, the array of time variant coefficients, based on \code{tParArray} using a coeffficients function. If no coefFunction provided, \code{ma1CoefArray} vector equals to \code{tParArray} vector.
#'
#' @param tParArray a T parameter, a.k.a. a sleeper parameter, a sequence of numbers {0, 1}.
#' @param coefFunction a coeffficients function, which is an underlying function used to generate \code{ma1CoefArray} vector.
#'
#' @return The vector of the array of time variant coefficients values corresponding to T parameter values.
#'
#' @examples
#'
#' createMA1CoefArray(coefFunction = sin, tParArray =  seq(from = 0, to = 1, length.out = 10))
#' createMA1CoefArray(coefFunction = sin, tParArray =  seq(from = 0, to = 1, by = 0.1))
#'
#'@export

createTVMA1CoefArray <- function(
                                 sampleSize
                                 )
{
  ma1CoefArray <- array(0, dim = sampleSize)

  for (sampleIndex in 1 : sampleSize)
  {

    ma1CoefArray [sampleIndex] <- customCoefFunction( sampleIndex / sampleSize)
  }
  return(ma1CoefArray)
}
