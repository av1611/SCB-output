source("R/normalDifferenceKernel.R")

  u = seq(-10, 10, 0.1)
  myK = normalDifferenceKernel (u)
  plot (myK ~ u)
  cat ("after plots")

