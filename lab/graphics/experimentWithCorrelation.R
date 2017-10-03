rm (list = ls())

lag <- 1
sampleSize <- 200


if(!dir.exists("../lab/graphics/out"))
{
  dir.create("../lab/graphics/out")
}

jpeg("../lab/graphics/out/mygraphWithCor.jpg")

myMA1CoefArray <- createTVMA1CoefArray(coefFunction = sin, sampleSize = sampleSize)

correlation <- computeCor(lag = lag, tvMa1CoefArray = myMA1CoefArray)

myTParArray <- createTParArray(tParCount = sampleSize)
plot(x = myTParArray, y = correlation, type="l")
dev.off()
