function1 = function (par1)
{
  cat("\n fileName=", par1)
}

mainFunction = function ()
{
  mySampleSize=10
  myLag = 1
  myBandwidth = 0.5
  myNonCoverageProbability=0.05
  fileName <- paste("ss", mySampleSize, "l", myLag, "bandW", myBandwidth, "alpha"
                    , myNonCoverageProbability, sep = "_")
  function1(fileName)
}
