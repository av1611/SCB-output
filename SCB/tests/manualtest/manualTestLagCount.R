maxGraph <- 1000

step <- 100

sampleSizes <- seq(1,maxGraph,by=step)
lags <- seq(1,10,by = 1)

lagCounts <- matrix(nrow = length(sampleSizes),ncol = length(lags))

for(i in 1:length(sampleSizes))
{
  for(j in 1:length(sampleSizes))
  {
    lagCounts[i,j]<- computeLagCount(sampleSize = sampleSizes[i],lag = lags[j])
  }
}

fileName <- "LagCount"

df <- data.frame(lagCounts)

saveCVS(fileName = fileName, doPath(), df)

saveJpg(fileName = fileName, doPath())

wireframe(lagCounts,row.values = sampleSizes,column.values = lags,
          xlab="sampleSizes",ylab="Lags")

dev.off()
