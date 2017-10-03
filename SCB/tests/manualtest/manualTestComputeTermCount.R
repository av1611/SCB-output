maxGraph <- 1000

step <- 1

sampleSizes <- seq(0,maxGraph,by=step)

termCounts <-numeric(length = length(sampleSizes))

for(i in 1:length(sampleSizes))
{
  termCounts[i] <- computeTermCount(i)
}

fileName <- "TermCount_SampleSize"

df <- data.frame(sampleSizes,termCounts)

saveJpg(fileName = fileName, doPath())

saveCVS(fileName = fileName, doPath(), df)

plot(sampleSizes,termCounts,type = "l")

title("TermCount and sampleSize")
dev.off()
