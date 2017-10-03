normalDifferenceKernelTestFunction = function () {
  cat("\n Testing \'normalDifferenceKernel\' properties \n")
  u = seq(-10, 10, 0.1)
  myK = normalDifferenceKernel (u)
  plot (myK ~ u)
}

# recollect the requirements to kernel
# like integration to 1
test_that("Testing normalDifferenceKernel properties",
          normalDifferenceKernelTestFunction())




