
rm (list = ls())
myPath <- "../lab/graphics/out"
if(!dir.exists(myPath))
{
  dir.create(myPath)
}
myFileName <- gsub(" ","_",paste("myGraphWithBand&Cor_",Sys.time(),".jpg",sep = ""))
myFileName <- gsub(":","_",myFileName)
jpeg(paste(myPath,"/",myFileName,sep=""))


mySampleSize=10
myTParCount = 10
mockTParArray <- createTParArray(tParCount = myTParCount)
mockTVMA1Array <- createTVMA1CoefArray(coefFunction = sin,sampleSize = mySampleSize)
myLag = 1
myLagCount = 4
myKernel = normalDifferenceKernel
myBandwidth = 0.5
myNonCoverageProbability=0.05

correlation <- computeCor(lag = myLag, coefFunction = sin, tParArray = mockTParArray)

band1 <- createBand(
  tParArray = mockTParArray,
  lag = myLag,
  lagCount = myLagCount,
  bandwidth = myBandwidth,
  kernel = myKernel,
  sampleSize = mySampleSize,
  nonCoverageProbability = myNonCoverageProbability )

minBand <- min(band1)
maxBand <- max(band1)
yMin = -1
yMax = 1
while(T)
{
  if(yMin < minBand)
  {
    if(yMax > maxBand)
    {
      break();
    }
    else
    {
      yMax = yMax + 0.01
    }
  }
  else
  {
    yMin = yMin - 0.01
  }
}
middle = (band1[,1] + band1[,2])/2

saveData <- data.frame(band1,middle,correlation)
myFileName <- gsub(" ","_",paste("myGraphWithBand&Cor_",Sys.time(),".csv",sep = ""))
myFileName <- gsub(":","_",myFileName)
write.csv(saveData, paste(myPath,"/",myFileName,sep=""))

plot(x=c(yMax:yMin),y=c(yMax:yMin),type = "n", xlim=c(0:1),main= "Band and Correlation")
lines(x=mockTParArray,y = band1[,1], type = "l", col="green")
lines(x=mockTParArray,y = band1[,2], type = "l", col="blue")
lines(x=mockTParArray,y = middle, type = "l", col = "brown")
lines(x=mockTParArray,y = correlation, type = "l", col="red")
lineArray <- c("Upper", "Middle", "Lower", "Correlation")
legend("bottom", title ="Correlation", lineArray, fill =c("blue", "brown", "green", "red"))
dev.off()
