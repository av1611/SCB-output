
maxGraph <- 1000

step <- 1

sampleSizes <- seq(0,maxGraph*10,by=step*10)
lags <- seq(0,maxGraph,by = step)

compatibles <- matrix(nrow = length(sampleSizes),ncol = length(lags))
cocoses <- numeric(length = length(sampleSizes))

for(i in 1:length(lags))
{
  for(j in 1:length(sampleSizes))
  {
    compatibles[i,j] <- isInputCompatible(sampleSize = sampleSizes[j], lag = lags[i])

    if(isInputCompatible(sampleSize = sampleSizes[j],lag = lags[i]))
    {
      cocoses[i] <- j
      break
    }
  }
}

fileName <- "IsinPutCompatible"

xMaxLim <- floor(maxGraph * 0.8)
yMaxLim <- max(cocoses)

saveJpg(fileName = fileName, doPath())
saveCVS(fileName = fileName, doPath(), compatibles)

plot(cocoses ~ lags,type = "l",xlim=c(1,xMaxLim),ylim=c(1,10000),
     col="red",ylab = "sampleSizes")
title("Is input Compatible")
#legend("bottomright",c("compatibles"),fill = c("red"))
#lines(cocoses,lags,col="red")

dev.off()


