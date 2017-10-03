#'
#'@export

saveCVS <- function(fileName, path, dataToSave)
{
  if(!dir.exists(path))
  {
    dir.create(path)
  }

  curTime <- Sys.time()
  curTime <- gsub(":", "-", curTime)
  curTime <- gsub(" ", "_", curTime)
  fileName <- gsub(" ","_", paste0(fileName, "_", curTime,".csv"))
  fileName <- gsub(":","_", fileName)
  write.csv(dataToSave, paste0(path, "/", fileName))
}
