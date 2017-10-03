
if (! require ("testthat"))
  install.packages ("testthat")
library ("testthat")

# make it work out of the box

test_file("tests/testthat/tests_01_createTParArray.R",
          reporter = "rstudio")
