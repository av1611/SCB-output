computeAllCorHatsTestFunction = function () {
  cat("\n Testing \'tests_33_computeAllCorHats\' \n")
  mySampleSize <- 62
  myTParCount = 10
  mockTParArray <- createTParArray(tParCount = myTParCount)
  # may be different
  mockTVMA1Array <- createTVMA1CoefArray(sampleSize = mySampleSize)

  mockSample <- createSample(sampleSize = mySampleSize)

  myLag = 3
  myLagCount <- computeLagCount(mySampleSize,myLag)
  myKernel = normalDifferenceKernel

  errorIfNotInputCompatible(alpha = 0.1,lag = myLag, sampleSize = mySampleSize)
  # form all rho hats
  # first try fakes
  Start=Sys.time()
  mockAllCorHats <- computeAllCorHats(tParArray = mockTParArray,
                                      lagCount = myLagCount,
                                      sample = mockSample)
  End=Sys.time()
  Duration=End-Start
cat("\nAllCorHats",mockAllCorHats)
cat("\nDuration",Duration,"\n")
cat("=====================")
cat("\nTest parameters :","\n")
cat("SampleSize= ",mySampleSize,"\n")
cat("TParCount= ",myTParCount,"\n")
cat("Lag= ",myLag,"\n")
cat("LagCount= ",myLagCount,"\n")
  # # print the first 9 elements. Hope, sizes are greater than 3.
  # lag = 1 # actual or real lag = 0
  # # mockCorHatsSecCol <- eval(parse(text = paste0("mockAllCorHats$`", lag, "`")))
  # mockCorHatsFirstCol <- eval(parse(text = paste0("mockAllCorHats[ , ", lag, "]")))
  # cat("The first col of mockCorHats = ", mockCorHatsFirstCol, "\n")
  #
  # lag = 2 # actual or real lag = 1
  # # mockCorHatsSecCol <- eval(parse(text = paste0("mockAllCorHats$`", lag, "`")))
  # mockCorHatsSecCol <- eval(parse(text = paste0("mockAllCorHats[ , ", lag, "]")))
  # cat("The second col of mockCorHats = ", mockCorHatsSecCol, "\n")


  # same but by steps does not work
  # colRefName = paste0("mockCorHats$`", lag, "`")
  # parsedColRefName = parse (text = colRefName)
  # evaluatedColumn = eval (parsedColRefName)
  #
  # cat("colRefName = ", colRefName, "\n")
  # cat("parsedColRefName = ", parsedColRefName, "\n")
  # cat("evaluatedColumn = ", evaluatedColumn, "\n")

}

test_that("Testing \'computeAllCorHats\'",{
          computeAllCorHatsTestFunction()
            cat(" End of test computeAllCorHats","\n")
            cat("=====================")
          })
# lagIndexToRefer <- eval(parse(text = paste0("allCorHats$`", lagIndex + 1, "`")))


