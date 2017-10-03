#' @export

doAll <- function(scbParams) {
  scb <- list()
  scb$tParArray <- createTParArray(tParCount = scbParams$tParCount)
  scb$tvma1CoefArray <- createTVMA1CoefArray(sampleSize = scbParams$sampleSize)
  scb$noise <- createNoise(sampleSize = scbParams$sampleSize,
                           mean = scbParams$mean,
                           sd = scbParams$sd)
  scb$sample <- createSample(sampleSize = scbParams$sampleSize)
  scb$k <- normalDifferenceKernel(u = scbParams$u)
  scb$bwidth <- computeB(sampleSize = scbParams$sampleSize)
  scb$multiplierCovarianceByKernel <-  createMultiplierCovarianceByKernel(kernel = scbParams$kernel,
                                                                          bandwidth = scb$bwidth,
                                                                          sampleSize = scbParams$sampleSize    )
  scb$bootstrapMultiplier <- createBootstrapMultiplier(kernel = scbParams$kernel,
                                                       bandwidth = scb$bwidth,
                                                       sampleSize = scbParams$sampleSize)

  scb$lagCount <-  computeLagCount(sampleSize = scbParams$sampleSize,
                                   lag = scbParams$lag)
  scb$covHat <- computeCovHat(tParArray = scb$tParArray,
                              lag = scbParams$lag,
                              sample = scb$sample)
  scb$corHat <- computeCorHat(tParArray = scb$tParArray,
                              lag = scbParams$lag,
                              sample = scb$sample)
  scb$allCorHats <- computeAllCorHats(tParArray = scb$tParArray,
                                      lagCount = scb$lagCount,
                                      sample = scb$sample)
  scb$betaLRVHat <- computeBetaLRVHat(tParArray = scb$tParArray,
                                      lag = scbParams$lag,
                                      sample = scb$sample,
                                      allCorHats = scb$allCorHats)
  scb$corArray <- computeCor(lag = scbParams$lag,
                             tParArray = scb$tParArray)
  scb$me <- computeMEbyCovHat(tParArray = scb$tParArray,
                              lag = scbParams$lag,
                              lagCount = scb$lagCount,
                              sample = scb$sample,
                              bandwidth = scb$bwidth,
                              nonCoverageProbability = scbParams$nonCoverageProbability,
                              allCorHats = scb$allCorHats,
                              C_K = scbParams$c_k,
                              PHI_K_NORMAL_DIFF = scbParams$phi_k_norm_diff)
  scb$band <- createBand(scb$tParArray,
                         scbParams$lag,
                         scb$lagCount,
                         scb$bwidth,
                         scbParams$sampleSize)
  scb$bandsBrick <- createBandsBrick(tParArray = scb$tParArray,
                                     lag = scbParams$lag,
                                     lagCount = scb$lagCount,
                                     bandwidth = scb$bwidth,
                                     kernel = normalDifferenceKernel,
                                     sampleSize = scbParams$sampleSize,
                                     nonCoverageProbability = scbParams$nonCoverageProbability,
                                     replicationCount = scbParams$replicationCount)
  scb$isCovered <- computeIsCovered(scb$band,
                                    scb$corArray,
                                    sampleSize = scbParams$sampleSize,
                                    bandwidth = scb$bwidth,
                                    lag = scbParams$lag,
                                    replicationCount = scbParams$replicationCount,
                                    superReplicationCount = scbParams$superReplicationCount)
  scb$isCoveredArray <- computeIsCoveredArray(bandsBrick = scb$bandsBrick,
                                              corArray = scb$corArray,
                                              sampleSize = scbParams$sampleSize,
                                              bandwidth = scb$bwidth,
                                              lag = scbParams$lag,
                                              replicationCount = scbParams$replicationCount,
                                              superReplicationCount = scbParams$SuperReplicationCount)
  scb$nonCoverageFreq <- computeNonCoverageFreq(replicationCount = scbParams$replicationCount,
                                                sampleSize = scbParams$sampleSize,
                                                lagCount = scb$lagCount,
                                                lag = scbParams$lag,
                                                tParArray = scb$tParArray,
                                                corArray = scb$corArray,
                                                kernel = scbParams$kernel,
                                                bandwidth = scb$bwidth,
                                                nonCoverageProbability = scbParams$nonCoverageProbability,
                                                superReplicationCount = scbParams$superReplicationCount)
  # scb$nonCoverageFreqArray <- computeNonCoverageFreqArray(superReplicationCount = scbParams$superReplicationCount,
  #                                                         replicationCount = scbParams$replicationCount,
  #                                                         sampleSize = scbParams$sample,
  #                                                         lag = scbParams$lag,
  #                                                         lagCount = scb$lagCount,
  #                                                         tParArray = scb$tParArray,
  #                                                         kernel = scbParams$kernel,
  #                                                         bandwidth = scb$bwidth,
  #                                                         nonCoverageProbability = scbParams$nonCoverageProbability)
  # scb$doubleAlphaArray = createDoubleAlphaArray(superReplicationCount = scbParams$superReplicationCount,
  #                                               replicationCount = scbParams$replicationCount,
  #                                               sampleSize = scbParams$sampleSize,
  #                                               # alphaArray =alphaArray,
  #                                               lag = scbParams$lag,
  #                                               lagCount = scb$lagCount,
  #                                               tParArray = scb$tParArray,
  #                                               kernel = scbParams$kernel,
  #                                               bandwidth = scb$bwidth)
  # saveBand(band = band,
  #          corArray = corArray,
  #          sampleSize = scbParams$sampleSize,
  #          lag = scbParams$lag,
  #          replicationCount = scbParams$replicationCount,
  #          bandwidth = scb$bwidth,
  #          superReplicationCount = scbParams$superReplicationCount)
  #
  # saveNonCoverageFreqArray(nonCoverageProbability = nonCoverageProbability,
  #                          alphaHatArray = scb$isCoveredArray,
  #                          sampleSize = scbParams$sampleSize,
  #                          replicationCount = scbParams$replicationCount,
  #                          bandwidth = scb$bwidth,
  #                          lag = scbParams$lag,
  #                          superReplicationCount = scbParams$superReplicationCount)
  # saveDoubleAplhaHatArray(nonCoverageProbabilities = scbParams$nonCoverageProbabilities,
  #                         alphaHats = alphaHats,
  #                         sampleSize = scbParams$sampleSize,
  #                         lag = scbParams$lag,
  #                         replicationCount = scbParams$replicationCount,
  #                         superReplicationCount = scbParams$superReplicationCount,
  #                         bandwidth = scb$bwidth)

  scb
}
