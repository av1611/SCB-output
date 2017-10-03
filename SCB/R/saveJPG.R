#'
#'@export

saveJpg <- function(fileName, path)
{
  if(!dir.exists(path))
  {
    dir.create(path)
  }

  fileName <- gsub("\\.","", fileName)
  curTime <- Sys.time()
  curTime <- gsub(":", "-", curTime)
  curTime <- gsub(" ", "_", curTime)
  fileName <- gsub(" ","_", paste0(fileName, "_", curTime,".jpg"))
  fileName <- gsub(":", "_", fileName)
  jpeg(paste0(path, "/", fileName))
}
