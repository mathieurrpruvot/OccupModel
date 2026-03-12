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
## 3. BUILD ADJACENCY MATRIX
############################################################

build_adjacency <- function(nrow, ncol){

  # Total number of grid cells
N <- nrow*ncol
  # Create empty NxN adjacency matrix
  # W[i,j] = 1 if j is a neighbour of i
W <- matrix(0,N,N)

index <- function(r,c) (r-1)*ncol + c

	for(r in 1:nrow){
		for(c in 1:ncol){

			i <- index(r,c)

			neighbors <- list(
			c(r-1,c),
			c(r+1,c),
			c(r,c-1),
			c(r,c+1),
			c(r-1,c-1),
			c(r+1,c-1),
			c(r+1,c+1),
			c(r-1,c+1)
			)

			for(nb in neighbors){

				rr <- nb[1]
				cc <- nb[2]

				if(rr>=1 & rr<=nrow & cc>=1 & cc<=ncol){

					j <- index(rr,cc)
					W[i,j] <- 1

				}
			}
		}
	}

return(W)

}

W <- build_adjacency(nrow_grid,ncol_grid)

############################################################
## BUILD NEIGHBOUR LISTS (FASTER THAN MATRIX MULTIPLICATION)
############################################################

build_neighbors <- function(W){
  
  N <- nrow(W)
  maxNbr <- max(rowSums(W))
  
  nbr <- matrix(1, N, maxNbr)
  numNbr <- rowSums(W)
  
  for(i in 1:N){
    
    ids <- which(W[i,]==1)
    
    nbr[i,] <- i
    nbr[i,1:length(ids)] <- ids
    
  }
  
  list(
    nbr=nbr,
    numNbr=numNbr,
    maxNbr=maxNbr
  )
}
nbr_info <- build_neighbors(W)

nbr <- nbr_info$nbr
numNbr <- nbr_info$numNbr
maxNbr <- nbr_info$maxNbr

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
alphaNbr = 0.1,
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

	nOccNbr <- W %*% z[,t-1]

	gamma <- plogis(true$alpha0 +
						true$alphaNbr*(nOccNbr/numNbr) +
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
alphaNbr ~ dnorm(0,0.001)
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

	  for(k in 1:maxNbr){
	    zNbr[i,t-1,k] <- z[nbr[i,k],t-1]
	  }
	  
	  nOccNbr[i,t-1] <- sum(zNbr[i,t-1,1:maxNbr])

		logit(gamma[i,t]) <- alpha0 +
		alphaNbr*(nOccNbr[i,t-1]/numNbr[i]) +
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
nbr=nbr,
numNbr=numNbr,
maxNbr=maxNbr
)

data <- c(detec,land)

z_inits <- matrix(1,N,T)

inits <- list(
  z=z_inits,
  beta0=0,beta1=0,beta2=0,
  alpha0=0,alphaNbr=0,alpha2=0,alpha3=0,alphaControl=0,
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
 "alpha0","alphaNbr","alpha2","alpha3","alphaControl",
 "delta0", "delta1",
 "theta0", "theta_adequate","theta_good",
 "eta0","eta_adequate","eta_good",
 "phi0","phi_old","phi_fresh",
 "z"
))

mcmc <- buildMCMC(conf)

cmcmc <- compileNimble(mcmc,project=model)

samples <- runMCMC(cmcmc,
niter=2000,
nburnin=1000,
nchains=1,
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

results_A_effN <- list()

for(sim in 1:nSim){

cat("Simulation",sim,"\n")

land <- simulate_landscape(forest_prop,smooth_strength)

z <- simulate_occupancy(land,true)

detec <- simulate_detection(z,prop_cell_scout_gap,prop_scout_replicate,
								prop_cell_drone_gap,prop_drone_replicate,
								pWeekActive,siteActiveProb,seasonGapProb,maxGapLength,siteWideWeekMiss,
								ViewMaps)

samples <- fit_model(land,detec)
 
    
    
    
results_A_effN[[sim]] <- samples

}


############################################################
## 11. PARAMETER RECOVERY
############################################################

mcmc_to_matrix <- function(samples){
  
  if(inherits(samples,"mcmc.list")){
    samp <- do.call(rbind, samples)
  } else {
    samp <- as.matrix(samples)
  }
  
  return(samp)
}

# Etract Posterior summaries

extract_results <- function(results, paramNames){
  
  nSim <- length(results)
  
  postMean <- matrix(NA,nSim,length(paramNames))
  lowerCI <- matrix(NA,nSim,length(paramNames))
  upperCI <- matrix(NA,nSim,length(paramNames))
  
  colnames(postMean) <- paramNames
  colnames(lowerCI) <- paramNames
  colnames(upperCI) <- paramNames
  
  for(s in 1:nSim){
    
    samp <- mcmc_to_matrix(results[[s]])
    
    for(p in paramNames){
      
      postMean[s,p] <- mean(samp[,p])
      lowerCI[s,p] <- quantile(samp[,p],0.025)
      upperCI[s,p] <- quantile(samp[,p],0.975)
      
    }
  }
  
  list(
    postMean=postMean,
    lowerCI=lowerCI,
    upperCI=upperCI
  )
  
}


# obtain bias, RMSE, rel RMSE, coverage

evaluate_performance <- function(summaryRes,true){
  
  postMean <- summaryRes$postMean
  lowerCI <- summaryRes$lowerCI
  upperCI <- summaryRes$upperCI
  
  paramNames <- colnames(postMean)
  trueVec <- unlist(true)
  
  nSim <- nrow(postMean)
  
  bias <- colMeans(postMean) - trueVec
  
  rmse <- sqrt(
    colMeans(
      (postMean - matrix(trueVec,nSim,length(trueVec),byrow=TRUE))^2
    )
  )
  
  relBias <- bias / abs(trueVec)
  relRMSE <- rmse / abs(trueVec)
  
  coverage <- numeric(length(paramNames))
  
  for(p in 1:length(paramNames)){
    
    coverage[p] <- mean(
      (trueVec[p] >= lowerCI[,p]) &
        (trueVec[p] <= upperCI[,p])
    )
    
  }
  
  data.frame(
    parameter=paramNames,
    true=trueVec,
    bias=bias,
    RMSE=rmse,
    relBias=relBias,
    relRMSE=relRMSE,
    coverage=coverage
  )
  
}

# run eval

paramNames <- names(true)

summaryRes_A_effN <- extract_results(results_A_effN,paramNames)

performance_A_effN <- evaluate_performance(summaryRes_A_effN,true)

print(performance_A_effN)


#parameter recovery plot
plot_recovery <- function(postMean,true){
  
  paramNames <- colnames(postMean)
  
  df <- data.frame(
    true=rep(unlist(true),each=nrow(postMean)),
    estimate=as.vector(postMean),
    parameter=rep(paramNames,each=nrow(postMean))
  )
  
  ggplot(df,aes(true,estimate))+
    geom_point(alpha=0.6)+
    geom_abline(slope=1,intercept=0,color="red")+
    facet_wrap(~parameter,scales="free")+
    theme_minimal()+
    labs(
      title="Parameter recovery",
      x="True value",
      y="Posterior mean estimate"
    )
  
}


plot_recovery(summaryRes_A_effN$postMean,true)

plot_coverage <- function(performance){
  
  ggplot(performance,
         aes(parameter,coverage))+
    geom_bar(stat="identity")+
    geom_hline(yintercept=0.95,
               linetype="dashed",
               color="red")+
    ylim(0,1)+
    theme_minimal()+
    labs(
      title="Coverage probability",
      y="Coverage"
    )
  
}

plot_coverage(performance_A_effN)

# Record the end time
end_time <- Sys.time()

# Calculate the elapsed time
time_taken_A_effN <- end_time - start_time

# Print the result (default unit is seconds)
print(time_taken_A_effN)
