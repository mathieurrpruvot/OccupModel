# Model A: Adjust final code for simulation study
# Model B: Adjust simulation study for distance kernel
# Model C: adjust B to add false positive
# model D: adjust C to add a sign presence process (zt=1 create signs, and z<t=1 create persistent signs)
# model E: adjust D to estimate placement

# model F: adjust D to esimtate control locations


#Compare the 3 versions of model A before finalizing into a consolidated version
# here, we will compare efficiency and coneistency of results for the 3 version over 3 simulations of 2000 iterations

source("ModelA_baseline.r")
source("ModelA_efficientLoop.r")
source("ModelA_efficientNeighbour.r")