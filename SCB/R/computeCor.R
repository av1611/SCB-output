#' @title computeCor
#'
#' @name computeCor
#'
#' @aliases computeRho
#'
#' @description Computing correlation: \code{cor} (\eqn{\rho}) given \code{tvMa1CoefArray}.
#'
#' @param lag the lag
#' @param tvMa1CoefArray the tvMa1CoefArray
#'
#' @return Array \code{cor_lag(i)}, counter \code{i} corresponds to \code{tvMa1CoefArray_i = 0, 1, ..., length(tvMa1CoefArray)}. length(tvMa1CoefArray) = sampleSize
#'
#' @examples
#' tParCount <- 10
#' lag <- 0
#' tParArray <- createTParArray(tParCount = tParCount)
#' tvMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, tParArray = tParArray)
#'
#'@export

computeCor <- function(lag,tParArray) {
  tParCount=length(tParArray)
  corArray <- array(0, dim = tParCount)

  if (lag == 1) {
    for(tParIndex in 1: tParCount){
    corArray[tParIndex] <-
      customCoefFunction(tParArray[tParIndex]) / (1 + customCoefFunction(tParArray[tParIndex]) ^ 2)
    }
  }

  if (lag == 0) {
    corArray <- rep(1,tParCount)
  }

  corArray
}
