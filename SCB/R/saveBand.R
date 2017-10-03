#'
#'@export

saveBand <- function(band,
                     corArray,
                     sampleSize,
                     replicationCount,
                     lag,
                     superReplicationCount,
                     nonCoverageProbability,
                     fileName="")
{
  bandwidth=computeB(sampleSize=sampleSize)
  #cat ("\n saveBand is called")
  # forming arrays to draw
  middle  <- (band[ , 1] + band[ , 2]) / 2
  consolidatedBands <- data.frame(lower  = band[, 1],
                                  middle = middle,
                                  upper  = band[, 2],
                                  corArray = corArray)
  globMin <- min(apply(consolidatedBands, 2, min))
  globMax <- max(apply(consolidatedBands, 2, max))
  globSd <- max(apply(consolidatedBands, 2, sd))
  yMinMargin <- globMin - 2 * globSd
  yMaxMargin <- globMax + 4 * globSd
  # forming filenames
  fileName <- paste(fileName,"BandAndCorrelation",sep = "_")
  fileName <- paste(fileName,
                    "ss", sampleSize,
                    "l", lag,
                    "bandW", round(bandwidth, digits = 1 ), sep = "_")
  fileName <- paste(fileName,
                    "alpha",nonCoverageProbability,
                    "repC", replicationCount,
                    "SrepC", superReplicationCount, sep = "_")

  xLab <- "tParArray"

  subTitle <- paste(xLab,
                    "\n sampleSize = ", sampleSize,
                    ", lag = ", lag,
                    ", bandwidth = ", round(bandwidth, digits = 1),
                    ",\n replicationCount = ", replicationCount,
                    ", superRep = ", superReplicationCount,
                    ", Alpha = ", nonCoverageProbability,
                    sep = "")
  xMinMargin <- 0
  xMaxMargin <- 1

  directory <- "bands"

  path <- doPath(directory = directory)

  # saving data to CSV
  consolidatedBands <- data.frame(consolidatedBands,corArray,sampleSize,bandwidth,
                                  replicationCount,
                                  lag,
                                  superReplicationCount,
                                  nonCoverageProbability)
  saveCVS(fileName = fileName, path = path, dataToSave = consolidatedBands)

  tParCount <- length(corArray)
  tParArray <- createTParArray(tParCount - 1)

  # opening a graphical device and saving the drawing to it
  saveJpg(fileName = fileName, path = path)
  plot(x = tParArray,
       y = consolidatedBands$lower,
       type = "l",
       col = "blue",
       xlim = c(xMinMargin, xMaxMargin),
       ylim = c(yMinMargin, yMaxMargin),
       xlab = "",
       ylab = "Correlation",
       yaxs = "i")
  title(main = "Band and Correlation", sub = subTitle)
  lines(x = tParArray, y = consolidatedBands$middle, type = "l", col = "brown")
  lines(x = tParArray, y = consolidatedBands$upper, type = "l", col = "green")
  lines(x = tParArray, y = consolidatedBands$corArray, type = "l", col = "red")
  lineArray <- c("Upper", "Middle", "Lower", "Correlation")
  legend("topright",
         title = "Correlation",
         lineArray,
         fill = c("green", "brown", "blue", "red"))
  graphics.off()
}
