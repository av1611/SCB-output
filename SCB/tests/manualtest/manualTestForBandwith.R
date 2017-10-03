sampleSize <- 10 : 1000
bandwidth<-computeB(sampleSize)
path<-doPath()
fileName<-"bandwidth_for_SampleSize"
saveJpg(fileName = fileName,path = path)
plot(bandwidth~sampleSize,type="l")
title("Bandwidth For SampleSize")
graphics.off()
