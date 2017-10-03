createTParArrayFunction <- function (tParCount=10) {
  cat ("\n Testing \'tests_01_createTParArray\' for for length and being numeric \n")
  tParCount = tParCount
  start=Sys.time()
  mockTParArray <- createTParArray(tParCount = tParCount)
  end=Sys.time()
  duration=end-start
  cat ("mockTParArray[1:5] =", mockTParArray[1:5], "\n",
       "length(mockTParArray) =", length(mockTParArray), "\n")
  cat("Duration= ",duration,"\n")
  cat("=====================")
  cat("\nTest parameters :","\n")
  cat("TparCount= ",tParCount,"\n")

  # The actual testing
  expect_that(mockTParArray, is_a("numeric"))
}



test_that("Testing \'createTParArray\' for being an array", {
  mockTParArray <- createTParArrayFunction()

cat("End of test createTParArray","\n")
cat("=====================")
  }
)
