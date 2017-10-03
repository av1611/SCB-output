cat("before")
quit("my error message upon stopping")
quit(save = "my error message upon stopping")
cat("after")


stop("my error message ")

myFunction = function ()
{
  cat ("\n before stop, inside myFunction \n")
  # stop("my error message inside myFunction")
  cat ("\n after stop, still inside myFunction \n")
  # return (15)
  return
}
cat ("\n before function call \n")
returnValue = myFunction()
returnValue
cat ("\n after function call \n")


quit(save = "mySave.Rdata")

