computeMEbyCovHatFunction = function ()
{

cat("\n Testing \'tests_42_computeMEbyCovHatFunction\' \n")


  myNonCoverageProbability = 0.05
  mySampleSize = 300
  myTParCount = 10

  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different

  mockSample <- createSample(sampleSize = mySampleSize)

  myLag = 1


  errorIfNotInputCompatible(alpha = myNonCoverageProbability,sampleSize = mySampleSize,lag = myLag)

  myLagCount = computeLagCount(lag=myLag,sampleSize = mySampleSize)


  myKernel = normalDifferenceKernel

  # form all rho hats
  # first try fakes
  mockAllCorHats <- computeAllCorHats(tParArray = mockTParArray,
                                      lagCount = myLagCount,
                                      sample = mockSample
                                      )
  Start=Sys.time()
  me <- computeMEbyCovHat(tParArray = mockTParArray,
                          lag = myLag,
                          lagCount = myLagCount,
                          sample = mockSample,

                          nonCoverageProbability = myNonCoverageProbability,
                          allCorHats = mockAllCorHats,
                          C_K = -1.978325,
                          # int_sq_der = 0.306951,
                          PHI_K_NORMAL_DIFF = 0.4065)

  End=Sys.time()
  Duration=End-Start

  # expect_that(mockBand, is_a("matrix"))
  # expect_that(dim(mockBand)[1], equals(2))  # the number of rows
  # expect_that(dim(mockBand)[2], equals(10)) # the number of cols

  cat("\nme = ", me, "\n")
  cat("Duration= ",Duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("SampleSize= ",mySampleSize,"\n")
  cat("TParCount= ",myTParCount,"\n")
  cat("Lag= ",myLag,"\n")
  cat("NonCoverageProbability= ",myNonCoverageProbability,"\n")

}

test_that("Testing \'computeMEbyCovHat\'", {
  computeMEbyCovHatFunction()
  cat(" End of test computeMEbyCovHat","\n")
  cat("=====================")
}
)

