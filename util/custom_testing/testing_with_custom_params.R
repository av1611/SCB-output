getwd()

source("../util/custom_testing/doAll.R")

devtools::load_all(".")

scbParams <- list(tParCount = 10,
                  sampleSize = 10,
                  mean = 0,
                  sd = 1,
                  kernel = normalDifferenceKernel,
                  u = seq(-10, 10, 0.1),
                  lag = 2,
                  nonCoverageProbability = 0.05,
                  c_k = -1.978325,
                  phi_k_norm_diff = 0.4065,
                  alphaCount = 10,
                  replicationCount = 10,
                  superReplicationCount = 10,
                  bwidth = 1)

scbList <- doAll(scbParams = scbParams)



