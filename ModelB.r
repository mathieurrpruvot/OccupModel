############################################################
## DYNAMIC MULTI-METHOD OCCUPANCY SIMULATION STUDY
## Wild Pig Monitoring
##
## Author: Mathieu Pruvot
## Purpose:
## 1) Simulate landscape
## 2) Simulate occupancy dynamics
## 3) Simulate detection data
## 4) Fit Nimble model
## 5) Evaluate parameter recovery

##
############################################################
##          MODEL B
############################################################


############################################################
## 1. LOAD LIBRARIES
############################################################

library(nimble)
library(terra)
library(ggplot2)
library(coda)

#set.seed(123)




############################################################
## 2. GLOBAL SETTINGS
############################################################

# Record the start time
start_time <- Sys.time()



nrow_grid <- 20
ncol_grid <- 20

N <- nrow_grid * ncol_grid
T <- 8

nScout <- 4
nDrone <- 4
maxCam <- 3
nWeek  <- 6

nSim <- 3   # number of simulation replicates


############################################################
## 3. BUILD DISTANCE MATRIX
############################################################

build_distance_matrix <- function(nrow,ncol){
  
  coords <- expand.grid(
    row=1:nrow,
    col=1:ncol
  )
  
  N <- nrow(coords)
  
  D <- matrix(0,N,N)
  
  for(i in 1:N){
    for(j in 1:N){
      
      D[i,j] <- sqrt(
        (coords$row[i]-coords$row[j])^2 +
          (coords$col[i]-coords$col[j])^2
      )
      
    }
  }
  
  return(D)
  
}

D <- build_distance_matrix(nrow_grid,ncol_grid)

############################################################
## BUILD DISPERSAL KERNEL
############################################################

sigma <- 2

K <- exp(-D/sigma)
diag(K) <- 0


############################################################
## 4. SIMULATE LANDSCAPE
############################################################

simulate_landscape <- function(forest_prop,smooth_strength){

# patchy forest




# ---- Create random noise ----
m <- matrix(runif(N), nrow = nrow_grid, ncol = ncol_grid)

# ---- Smooth using focal averaging ----


rasterforest <- rast(m)

# Create smoothing kernel
smooth_kern <- matrix(1, nrow = smooth_strength, ncol = smooth_strength)

# Apply smoothing
r_smooth <- focal(rasterforest, w = smooth_kern, fun = mean, na.rm = TRUE)

# Convert to binary forest/non-forest
threshold <- quantile(values(r_smooth), probs = 1 - forest_prop)
forest <- r_smooth > threshold
forest <- c(as.matrix(forest)*1)

distCrop <- as.vector(scale(runif(N)))

control <- matrix(rbinom(N*T,1,0.7),N,T)

list(
forest=forest,
distCrop=distCrop,
control=control
)

}


############################################################
## 5. TRUE PARAMETERS
############################################################

true <- list(

beta0 = -1,
beta1 = 1.5,
beta2 = 0,

alpha0 = -2,
alphaSpread = 0.15,
alpha2 = 1,
alpha3 = 0,
alphaControl = -1.5,

delta0 = -1,
delta1 = 1,

theta0 = 0,
theta_adequate = 0.4,
theta_good = 0.8,

eta0 = 0,
eta_adequate = 0.5,
eta_good = 1,

phi0 = 0,
phi_old = 0.5,
phi_fresh = 1
)


############################################################
## 6. SIMULATE OCCUPANCY
############################################################

simulate_occupancy <- function(land,true){

forest <- land$forest
distCrop <- land$distCrop
control <- land$control

z <- matrix(0,N,T)

psi1 <- plogis(true$beta0 +
true$beta1*forest +
true$beta2*distCrop)

z[,1] <- rbinom(N,1,psi1)

for(t in 2:T){

  spread <- K %*% z[,t-1]
  
  gamma <- plogis(
            true$alpha0 +
            true$alphaSpread*spread +
						true$alpha2*forest +
						true$alpha3*distCrop +
						true$alphaControl*control[,t])

	epsilon <- plogis(true$delta0 +
						true$delta1*control[,t])

	psi <- z[,t-1]*(1-epsilon) +
			(1-z[,t-1])*gamma

	z[,t] <- rbinom(N,1,psi)

}
#plot the change in occupancy over time
dftrueocc <- data.frame(occ=colSums(z),time=c(1:T))
ggplot(dftrueocc) +
  geom_line(aes(time,occ))

return(z)


}


############################################################
## 7. SIMULATE DETECTION DATA
############################################################







simulate_detection <- function(z,prop_cell_scout_gap,prop_scout_replicate,
								prop_cell_drone_gap,prop_drone_replicate,
								pWeekActive,siteActiveProb,seasonGapProb,maxGapLength,siteWideWeekMiss,
								ViewMaps){

### detection covariates

scoutCond <- array(sample(1:3,N*T*nScout,TRUE),
					c(N,T,nScout))

droneCond <- array(sample(1:3,N*T*nDrone,TRUE),
					c(N,T,nDrone))

camFOV <- array(sample(1:3,N*T*maxCam*nWeek,TRUE),
					c(N,T,maxCam,nWeek))


### SCOUT DATA

monitorScout <- matrix(1, N, T)

for(i in 1:N){
  
  if(runif(1) < prop_cell_scout_gap){        # 30% of cells experience a monitoring gap
    
    start_gap <- sample(1:(T-3),1)
    length_gap <- sample(2:3,1)
    
    monitorScout[i, start_gap:(start_gap+length_gap-1)] <- 0
  }
}


yScout <- array(NA,c(N,T,nScout))

for(i in 1:N){
  for(t in 1:T){
    
    if(monitorScout[i,t] == 1){
      
      for(j in 1:nScout){
        
        if(runif(1) < prop_scout_replicate){  # replicate conducted
          cond <- scoutCond[i,t,j]
          p <- plogis(true$theta0 +
            (cond==2)*true$theta_adequate +
            (cond==3)*true$theta_good)
          yScout[i,t,j] <- rbinom(1,1,z[i,t]*p)
        }
        
      }
      
    }
    
  }
}



### DRONE DATA

#prop_cell_drone_gap,prop_drone_replicate

monitorDrone <- matrix(1, N, T)

for(i in 1:N){
  if(runif(1) < prop_cell_drone_gap){        # 30% of cells experience a monitoring gap
    
    start <- sample(1:(T-5),1)
    length_gap <- sample(2:5,1)
    
    monitorDrone[i, start:(start+length_gap-1)] <- 0
  }
}

yDrone <- array(NA, c(N,T,nDrone))
for(i in 1:N){
  for(t in 1:T){
    
    if(monitorDrone[i,t] == 1){
      
    
    for(f in 1:nDrone){
      
      if(runif(1) < prop_drone_replicate){
        cond <- droneCond[i,t,f]
        p <- plogis(true$eta0 +
          (cond==2)*true$eta_adequate +
          (cond==3)*true$eta_good)
        yDrone[i,t,f] <- rbinom(1,1,z[i,t]*p)
    }
    }
  }
  }
}


### CAMERA DATA

# Initialize arrays
yCam <- array(NA, c(N, T, maxCam, nWeek))  # final detection array
nCam <- matrix(0, N, T)                    # cameras deployed per site-season
periodActive <- matrix(1, N, T)            # whether a season is monitored
weekActive <- array(0, c(N, T, maxCam, nWeek))  # weekly replicate activity

# Step 1: Determine site-level monitoring
siteActive <- rbinom(N, 1, siteActiveProb)  # some sites never monitored

# Step 2: Generate seasonal gaps within active sites
for(i in 1:N){
  if(siteActive[i]==1 && runif(1)<seasonGapProb){
    start <- sample(1:(T-maxGapLength), 1)
    gapLength <- sample(2:maxGapLength,1)
    periodActive[i, start:(start+gapLength-1)] <- 0
  }
}

# Step 3: Simulate cameras and weekly replicates
for(i in 1:N){
  
  if(siteActive[i]==1){   # skip completely inactive sites
    
    for(t in 1:T){
      
      if(periodActive[i,t]==1){  # skip missing seasons
        
        # Variable number of cameras per site-season
        nCam[i,t] <- sample(1:maxCam, 1)
        
        for(k in 1:nCam[i,t]){
          
          for(w in 1:nWeek){
            
            # site-wide week missing
            if(runif(1) < siteWideWeekMiss){
              weekActive[i,t,k,w] <- 0
            } else {
              # camera-level weekly probability
              weekActive[i,t,k,w] <- rbinom(1,1,pWeekActive)
            }
            
            # Assign detection if active, otherwise leave NA
            if(weekActive[i,t,k,w]==1){
              fov <- camFOV[i,t,k,w]
              p <- plogis(true$phi0 +
                (fov==2)*true$phi_old +
                (fov==3)*true$phi_fresh)
              yCam[i,t,k,w] <- rbinom(1,1,z[i,t]*p)
            }
            
          } # end weeks
        } # end cameras
      } # end period active
    } # end seasons
  } # end site active
} # end sites

# Visualiztion of camera trapping
#####
if(ViewMaps){


effort <- matrix(0, nrow_grid, ncol_grid)

for(i in 1:N){
  
  r <- ((i-1) %/% ncol_grid) + 1
  c <- ((i-1) %% ncol_grid) + 1
  
  effort[r,c] <- sum(nCam[i,])
}

df <- expand.grid(
  row=1:nrow_grid,
  col=1:ncol_grid
)

df$effort <- as.vector(effort)

ggplot(df, aes(col,row,fill=effort)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title="Total Camera Effort per Grid Cell",
    x="Column",
    y="Row"
  ) +
  theme_minimal()

siteCoverage <- rowSums(nCam>0)

coverageMat <- matrix(siteCoverage, nrow_grid, ncol_grid)

df$coverage <- as.vector(coverageMat)

ggplot(df, aes(col,row,fill=coverage)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_equal() +
  scale_y_reverse() +
  labs(
    title="Number of Seasons Monitored per Cell"
  ) +
  theme_minimal()

weekEffort <- matrix(0, nrow_grid, ncol_grid)

for(i in 1:N){
  
  r <- ((i-1) %/% ncol_grid) + 1
  c <- ((i-1) %% ncol_grid) + 1
  
  weekEffort[r,c] <- sum(weekActive[i,1,,])
}

df$weekEffort <- as.vector(weekEffort)

ggplot(df, aes(col,row,fill=weekEffort)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_equal() +
  scale_y_reverse() +
  labs(title="Camera Weeks Active (Season 1)") +
  theme_minimal()
}



list(
scoutCond=scoutCond,
droneCond=droneCond,
camFOV=camFOV,
yScout=yScout,
yDrone=yDrone,
yCam=yCam

)

}


############################################################
## 8. NIMBLE MODEL
############################################################

code <- nimbleCode({

beta0 ~ dnorm(0,0.001)
beta1 ~ dnorm(0,0.001)
beta2 ~ dnorm(0,0.001)

alpha0 ~ dnorm(0,0.001)
alphaSpread ~ dnorm(0,0.001)
alpha2 ~ dnorm(0,0.001)
alpha3 ~ dnorm(0,0.001)
alphaControl ~ dnorm(0,0.001)

delta0 ~ dnorm(0,0.001)
delta1 ~ dnorm(0,0.001)

theta0 ~ dnorm(0,0.01)
theta_adequate ~ dnorm(0,0.01)
theta_good ~ dnorm(0,0.01)

eta0 ~ dnorm(0,0.01)
eta_adequate ~ dnorm(0,0.01)
eta_good ~ dnorm(0,0.01)

phi0 ~ dnorm(0,0.01)
phi_old ~ dnorm(0,0.01)
phi_fresh ~ dnorm(0,0.01)


for(i in 1:N){

	logit(psi[i,1]) <- beta0 +
		beta1*forest[i] +
		beta2*distCrop[i]

	z[i,1] ~ dbern(psi[i,1])


	for(t in 2:T){

	  spread[i,t-1] <- inprod(K[i,1:N],z[1:N,t-1])

		logit(gamma[i,t]) <- alpha0 +
		alphaSpread*spread[i,t-1] +
		alpha2*forest[i] +
		alpha3*distCrop[i] +
		alphaControl*control[i,t]

		logit(epsilon[i,t]) <- delta0 +
		delta1*control[i,t]

		psi[i,t] <- z[i,t-1]*(1-epsilon[i,t]) +
			(1-z[i,t-1])*gamma[i,t]

		z[i,t] ~ dbern(psi[i,t])

	}


	for(t in 1:T){

		for(j in 1:nScout){

			logit(pS[i,t,j]) <- theta0 +
			equals(scoutCond[i,t,j],2)*theta_adequate +
			equals(scoutCond[i,t,j],3)*theta_good

			yScout[i,t,j] ~ dbern(z[i,t]*pS[i,t,j])
		}


		for(f in 1:nDrone){

			logit(pD[i,t,f]) <- eta0 +
			equals(droneCond[i,t,f],2)*eta_adequate +
			equals(droneCond[i,t,f],3)*eta_good

			yDrone[i,t,f] ~ dbern(z[i,t]*pD[i,t,f])
		}


		for(k in 1:maxCam){
			for(w in 1:nWeek){

			logit(pC[i,t,k,w]) <- phi0 +
			equals(camFOV[i,t,k,w],2)*phi_old +
			equals(camFOV[i,t,k,w],3)*phi_fresh

			yCam[i,t,k,w] ~ dbern(z[i,t]*pC[i,t,k,w])

		}
		}

	}

}

})


############################################################
## 9. MODEL FIT FUNCTION
############################################################

fit_model <- function(land,detec){

constants <- list(
N=N,T=T,
nScout=nScout,
nDrone=nDrone,
maxCam=maxCam,
nWeek=nWeek,
K=K
)

data <- c(detec,land)

z_inits <- matrix(1,N,T)

inits <- list(
  z=z_inits,
  beta0=0,beta1=0,beta2=0,
  alpha0=0,alphaSpread=0,alpha2=0,alpha3=0,alphaControl=0,
  delta0=0,delta1=0,
  theta0=0,theta_adequate=0,theta_good=0,
  eta0=0,eta_adequate=0,eta_good=0,
  phi0=0,phi_old=0,phi_fresh=0
  )

model <- nimbleModel(code,
					constants=constants,
					data=data,
					inits=inits)

cmodel <- compileNimble(model)

conf <- configureMCMC(model, monitors=c(
  "beta0","beta1","beta2",
 "alpha0","alphaSpread","alpha2","alpha3","alphaControl",
 "delta0", "delta1",
 "theta0", "theta_adequate","theta_good",
 "eta0","eta_adequate","eta_good",
 "phi0","phi_old","phi_fresh",
 "z"
))

mcmc <- buildMCMC(conf)

cmcmc <- compileNimble(mcmc,project=model)

samples <- runMCMC(cmcmc,
niter=3000,
nburnin=1000,
nchains=2,
samplesAsCodaMCMC=TRUE)

return(samples)

}


############################################################
## 10. SIMULATION LOOP
############################################################

# Parameters

forest_prop <- 0.2      # target proportion forest
smooth_strength <- 7    # higher = larger patches
pWeekActive <- 0.85     # probability a weekly replicate is active
siteActiveProb <- 0.8   # probability a site is ever monitored
seasonGapProb <- 0.4    # probability of a seasonal gap at a site
maxGapLength <- 3       # max consecutive seasons missing
siteWideWeekMiss <- 0.05 # chance all cameras at site miss a week
prop_cell_scout_gap <- 0.6
prop_scout_replicate <- 0.7
prop_cell_drone_gap <- 0.6
prop_drone_replicate <- 0.75
ViewMaps <- F

results <- list()

for(sim in 1:nSim){

cat("Simulation",sim,"\n")

land <- simulate_landscape(forest_prop,smooth_strength)

z <- simulate_occupancy(land,true)

detec <- simulate_detection(z,prop_cell_scout_gap,prop_scout_replicate,
								prop_cell_drone_gap,prop_drone_replicate,
								pWeekActive,siteActiveProb,seasonGapProb,maxGapLength,siteWideWeekMiss,
								ViewMaps)

samples <- fit_model(land,detec)

  
post <- summary(samples)$statistics[,"Mean"]

results[[sim]] <- post

}


############################################################
## 11. PARAMETER RECOVERY
############################################################

paramNames <- names(true)

estimates <- matrix(NA,nSim,length(paramNames))

colnames(estimates) <- paramNames

for(i in 1:nSim){

for(p in paramNames){

estimates[i,p] <- results[[i]][p]

}
}

bias <- colMeans(estimates) - unlist(true)

rmse <- sqrt(colMeans((estimates -
matrix(unlist(true),
nSim,length(true),byrow=TRUE))^2))


performance <- data.frame(
parameter=paramNames,
true=unlist(true),
bias=bias,
RMSE=rmse
)

print(performance)


# Record the end time
end_time <- Sys.time()

# Calculate the elapsed time
time_taken <- end_time - start_time

# Print the result (default unit is seconds)
print(time_taken)