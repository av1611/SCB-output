#' @title createTParArray
#'
#' @name createTParArray
#'
#' @aliases createSleeper
#'
#' @description Create \code{T parameter array}.
#'
#' @param tParCount the desired length of \code{T parameter} vector.
#'
#' @return The vector of \code{T parameter array}.
#'
#' @examples
#' tParCount = 10
#' createTParArray(tParCount = tParCount)
#'
#'@export

createTParArray <- function(tParCount)
{
  tParArray <- seq(from = 0,
                   to   = 1 - 1 / tParCount,
                   by   = 1 / tParCount)
  tParArray <- c(tParArray,1)
}
