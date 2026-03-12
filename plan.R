# Model A: Adjust final code for simulation study
# Model B: Adjust simulation study for distance kernel
# Model C: adjust B to add false positive
# model D: adjust C to add a sign presence process (zt=1 create signs, and z<t=1 create persistent signs)
# model E: adjust D to estimate placement

# model F: adjust D to esimtate control locations


build_neighbors <- function(W){
  
  N <- nrow(W)
  
  maxNbr <- max(rowSums(W))
  
  nbr <- matrix(0,N,maxNbr)
  numNbr <- rowSums(W)
  
  for(i in 1:N){
    
    ids <- which(W[i,]==1)
    
    nbr[i,1:length(ids)] <- ids
    
  }
  
  list(nbr=nbr, numNbr=numNbr, maxNbr=maxNbr)
  
}