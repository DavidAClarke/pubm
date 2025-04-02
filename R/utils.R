## Loads an RData file, and returns it
loadRData <- function(fileName){
  
  load(fileName)
  
  get(ls()[ls() != "fileName"])
  
}
