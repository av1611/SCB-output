  createNonCoverageFreqDoubleArrayFunction <- function () {
    cat("\n Testing \'tests_65_createNonCoverageFreqDoubleArray\' \n")
    myTParCount  <- 10
    myTParArray  <- createTParArray(tParCount = myTParCount)
    mySuperReplicationCount <- 3
    myReplicationCount <- 6
    mySampleSize <- 100
    myLag <- 2
    myLagCount <- computeLagCount(lag = myLag,sampleSize = mySampleSize)


    nonCoverageProbabilities <- c(0.2,0.4,0.6,0.8)
    Start=Sys.time()
    doubleAlphaArray = createDoubleAlphaArray(
      superReplicationCount = mySuperReplicationCount,
      replicationCount = myReplicationCount,
      sampleSize = mySampleSize,
      alphaArray =nonCoverageProbabilities,
      lag = myLag,
      lagCount = myLagCount,
      tParArray = myTParArray,
      fileName = "tests_65_createNonCoverageFreqDoubleArray")
    End=Sys.time()
    duration=End-Start
    cat("\nDoubleAlphaHatArray= ",doubleAlphaArray)
    cat("\nAlphaArray: ",nonCoverageProbabilities)
    cat("\n size of double array= ",length(doubleAlphaArray))
    cat("\n Duration= ",duration,"\n")
    cat("=====================")
    cat("\nTest parameters :","\n")
    cat("SampleSize= ",mySampleSize,"\n")
    cat("TParCount= ",myTParCount,"\n")
    cat("Lag= ",myLag,"\n")
    cat("LagCount= ",myLagCount,"\n")
    cat("ReplicationCount= ",myReplicationCount,"\n")
    cat("SuperReplicationCount= ",mySuperReplicationCount,"\n")
    # expect_that(mockBand, is_a("matrix"))
    # expect_that(dim(mockBand)[1], equals(2))  # the number of rows
    # expect_that(dim(mockBand)[2], equals(10)) # the number of cols

  }

  test_that("Testing \' createNonCoverageFreqDoubleArray\'", {
    createNonCoverageFreqDoubleArrayFunction()
  })
