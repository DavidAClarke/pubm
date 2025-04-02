## Loads an RData file, and returns it
loadRData <- function(fileName){
  
  load(fileName)
  
  get(ls()[ls() != "fileName"])
  
}

## Simulate Poisson data within bounds
rlimpois <- function(n, lambda, lowlimit, toplimit){
  
  sample(x = lowlimit:toplimit, 
         size = n,
         prob = dpois(lowlimit:toplimit, lambda), 
         replace = TRUE)
  
}