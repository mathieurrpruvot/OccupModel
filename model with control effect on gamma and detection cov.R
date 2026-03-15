####
#Notes:
# This seems to work, but needs a few improvements:
#   - DONE: need to make detection proba more realistic (lower)
#   - DONE: need to make colonization dependent on control activities, to limit the boom-bust dynamics on longer time scales (DONE)
#   - need to switch over to the model allowing detection covariates
#   - need to start introducing data gaps to allow for non-sampling of certain cells





############################################################
## 1. SETUP
############################################################

library(nimble)
library(ggplot2)

set.seed(123)

nrow_grid <- 20
ncol_grid <- 20
N <- nrow_grid * ncol_grid
T <- 8

nScout <- 4
nDrone <- 4
nWeek  <- 6
maxCam <- 3

############################################################
## 2. BUILD ADJACENCY MATRIX
############################################################

build_adjacency <- function(nrow, ncol){
  
  # Total number of grid cells
  N <- nrow * ncol
  
  # Create empty NxN adjacency matrix
  # W[i,j] = 1 if j is a neighbour of i
  W <- matrix(0, N, N)
  
  # Function converting (row, column) to single index
  # Example: (r=2,c=3) in 10-column grid → cell number
  index <- function(r,c) (r-1)*ncol + c
  
  for(r in 1:nrow){
    for(c in 1:ncol){
      
      # Current cell index
      i <- index(r,c)
      
      # Define potential 8-neighbours
      neighbors <- list(
        c(r-1,c),  # up
        c(r+1,c),  # down
        c(r,c-1),  # left
        c(r,c+1),   # right
        c(r-1,c-1), #top left corner
        c(r+1,c-1), #bottom left corner
        c(r+1,c+1), #bottom right corner
        c(r-1,c+1)  #top right corner
      )
      
      for(nb in neighbors){
        rr <- nb[1]; cc <- nb[2]
        
        # Only keep valid cells (avoid edges going out of grid)
        if(rr>=1 & rr<=nrow & cc>=1 & cc<=ncol){
          
          j <- index(rr,cc)
          
          # Mark j as neighbour of i
          W[i,j] <- 1
        }
      }
    }
  }
  
  return(W)
}

W <- build_adjacency(nrow_grid, ncol_grid)

############################################################
## 3. SIMULATE ECOLOGICAL COVARIATES
############################################################

#forest    <- runif(N, 0, 1)

########################
# Create patchy forest #
########################


# ---- Parameters ----
set.seed(123)
library(terra)

n <- 20                 # grid size
forest_prop <- 0.2      # target proportion forest
smooth_strength <- 7    # higher = larger patches


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

# ---- Plot ----
plot(forest,
     col = c("lightgoldenrod", "darkgreen"),
     legend = FALSE,
     main = "Simulated Patchy Forest Habitat")

forest <- c(as.matrix(forest)*1)


distCrop  <- runif(N, 0, 1)
control   <- matrix(rbinom(N*T,1,0.7), N, T)

# 4. SIMULATE DETECTION COVARIATES

scoutCond <- array(sample(1:3, N*T*nScout, replace=TRUE),
                   c(N,T,nScout))

droneCond <- array(sample(1:3, N*T*nDrone, replace=TRUE),
                   c(N,T,nDrone))

camFOV <- array(sample(1:3, N*T*maxCam*nWeek, replace=TRUE),
                c(N,T,maxCam,nWeek))

############################################################
## 4. TRUE PARAMETERS
############################################################

true <- list(
  beta0 = -1,
  beta1 = 1.5,
  beta2 = 0,
  
  alpha0 = -2,
  alphaNbr = 0.1,
  alpha2 = 1,
  alpha3 = 0,
  alphaControl=-1.5,
  
  delta0 = -1,
  delta1 = 1,
  
  theta0 = 0, theta_adequate=0.4, theta_good=0.8,
  eta0   = 0,eta_adequate=0.5, eta_good=1,
  phi0   = 0, phi_old=0.5, phi_fresh=1
)

############################################################
## 5. SIMULATE LATENT OCCUPANCY
############################################################

z <- matrix(0, N, T)

psi1 <- plogis(true$beta0 +
                 true$beta1*forest +
                 true$beta2*distCrop)

z[,1] <- rbinom(N,1,psi1)

for(t in 2:T){
  
  nOccNbr <- W %*% z[,t-1]
  
  gamma <- plogis(true$alpha0 +
                    true$alphaNbr*nOccNbr +
                    true$alpha2*forest +
                    true$alpha3*distCrop +
                    true$alphaControl*control[,t])
  
  epsilon <- plogis(true$delta0 +
                      true$delta1*control[,t])
  
  psi <- z[,t-1]*(1-epsilon) +
    (1-z[,t-1])*gamma
  
  z[,t] <- rbinom(N,1,psi)
}

dftrueocc <- data.frame(occ=colSums(z),time=c(1:8))

ggplot(dftrueocc) +
  geom_line(aes(time,occ))

############################################################
## 6. SIMULATE DETECTION DATA
############################################################

## SCOUTING
#yScout <- array(0, c(N,T,nScout))
#for(i in 1:N){
#  for(t in 1:T){
#    for(j in 1:nScout){
#      p <- plogis(true$theta0)
#      yScout[i,t,j] <- rbinom(1,1,z[i,t]*p)
#    }
#  }
#}
#

monitorScout <- matrix(1, N, T)

for(i in 1:N){
  
  if(runif(1) < 0.3){        # 30% of cells experience a monitoring gap
    
    start <- sample(1:(T-3),1)
    length_gap <- sample(2:3,1)
    
    monitorScout[i, start:(start+length_gap-1)] <- 0
  }
}

yScout <- array(NA, c(N,T,nScout))

for(i in 1:N){
  for(t in 1:T){
    
    if(monitorScout[i,t] == 1){
      
      for(j in 1:nScout){
        
        if(runif(1) < 0.8){  # replicate conducted
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

## DRONE

monitorDrone <- matrix(1, N, T)

for(i in 1:N){
  if(runif(1) < 0.6){        # 30% of cells experience a monitoring gap
    
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
      
      if(runif(1) < 0.8){
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

## CAMERA
#nCam <- rep(maxCam, N)
#yCam <- array(0, c(N,T,maxCam,nWeek))
#
#for(i in 1:N){
#  for(t in 1:T){
#    for(k in 1:maxCam){
#      for(w in 1:nWeek){
#        p <- plogis(true$phi0)
#        yCam[i,t,k,w] <- rbinom(1,1,z[i,t]*p)
#      }
#    }
#  }
#}

############################################################
## CAMERA MONITORING WITH REALISTIC MISSINGNESS
############################################################


# Parameters
pWeekActive <- 0.85     # probability a weekly replicate is active
siteActiveProb <- 0.8   # probability a site is ever monitored
seasonGapProb <- 0.4    # probability of a seasonal gap at a site
maxGapLength <- 3       # max consecutive seasons missing
siteWideWeekMiss <- 0.05 # chance all cameras at site miss a week

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


####
# Visualiztion of camera trapping
#####

library(ggplot2)

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



############################################################
## 7. NIMBLE MODEL
############################################################

code <- nimbleCode({
  
  # Priors
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
      
      nOccNbr[i,t-1] <- inprod(W[i,1:N], z[1:N,t-1])
      
      logit(gamma[i,t]) <- alpha0 +
        alphaNbr*nOccNbr[i,t-1]/8 +
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
      
      for(k in 1:maxCam){                # only deployed cameras
        for(w in 1:nWeek){
          logit(pC[i,t,k,w]) <- phi0 +
            equals(camFOV[i,t,k,w],2)*phi_old +
            equals(camFOV[i,t,k,w],3)*phi_fresh
          yCam[i,t,k,w] ~ dbern(z[i,t] * pC[i,t,k,w])
          
        }
      }
      
    }
  }
})

############################################################
## 8. BUILD MODEL
############################################################

constants <- list(
  N=N, T=T,
  nScout=nScout,
  nDrone=nDrone,
  nWeek=nWeek,
  maxCam=maxCam,
  W=W
)

data <- list(
  yScout=yScout,
  yDrone=yDrone,
  yCam=yCam,
  forest=forest,
  distCrop=distCrop,
  control=control,
  scoutCond=scoutCond,
  droneCond=droneCond,
  camFOV=camFOV
)

z_init <- matrix(1,N,T)

inits <- list(
  z=z_init,
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
cmcmc <- compileNimble(mcmc, project=model)

samples <- runMCMC(cmcmc, niter=2000, nburnin=1000)

summary(samples[,c(1:19)])
traceplot(samples[,c(1:19)])

library(coda)
# Assuming 'samples' is a list of matrices from multiple chains
library(MCMCvis)
MCMCtrace(samples[,c(1:19)],open_pdf = F)
MCMCsummary(samples[,c(1:19)])
MCMCplot(samples)

#extract posterieor occupacy proba

library(coda)

# Extract z samples
z_samples <- samples[, grep("^z\\[", colnames(samples))]

colSums(matrix(colMeans(z_samples),nrow=N,ncol=8))

est_z <- matrix(colMeans(z_samples),nrow=N,ncol=8)

# Convert to matrix form
n_iter <- nrow(z_samples)

# Rebuild into array: iterations × N × T
z_array <- array(NA, c(n_iter, N, T))

counter <- 1
for(i in 1:N){
  for(t in 1:T){
    z_array[,i,t] <- z_samples[,counter]
    counter <- counter + 1
  }
}

# Posterior mean occupancy
psi_hat <- apply(z_array, c(2,3), mean)



#Map occupancy through time


to_grid <- function(x){
  matrix(x, nrow=nrow_grid, ncol=ncol_grid, byrow=TRUE)
}

library(ggplot2)
library(reshape2)
library(patchwork)

plot_grid <- function(mat, title){
  
  df <- melt(mat)
  colnames(df) <- c("row","col","value")
  
  ggplot(df, aes(col, row, fill=value)) +
    geom_tile() +
    scale_fill_viridis_c(option="C") +
    coord_fixed() +
    scale_y_reverse() +
    theme_minimal() +
    labs(title=title, fill="Occupancy")
}

plots <- list()

for(t in 1:T){
  
  true_map <- to_grid(z[,t])
  est_map  <- to_grid(est_z[,t])
  
  p1 <- plot_grid(true_map, paste("True Occupancy - t =", t))
  p2 <- plot_grid(est_map, paste("Posterior Mean - t =", t))
  
  plots[[t]] <- p1 + p2
}

wrap_plots(plots)


#Detection vs posterio fit

obs_detect <- apply(yScout, c(1,2), max)

fit_df <- data.frame(
  observed = as.vector(obs_detect),
  est_z  = as.vector(est_z)
)

ggplot(fit_df, aes(est_z, observed)) +
  geom_jitter(height=0.05, width=0) +
  geom_smooth(method="glm", method.args=list(family="binomial")) +
  theme_minimal() +
  labs(x="Posterior Occupancy",
       y="Observed Detection")

fit_df$observed == fit_df$est_z

#Invasion Wave visualization

true_occ <- colSums(z)
est_occ  <- colSums(matrix(colMeans(z_samples),nrow=N,ncol=8))

df_occ <- data.frame(
  time=1:T,
  True=true_occ,
  Estimated=est_occ
)

ggplot(df_occ) +
  geom_line(aes(time, True), size=1.2) +
  geom_line(aes(time, Estimated), linetype="dashed", size=1.2) +
  theme_minimal() +
  labs(y="Number of Occupied Cells",
       title="Invasion Dynamics: Simulated (line) vs. Estimated (dash)")


#Effect of number of occupied neighbours

alphaNbr_post <- samples[,"alphaNbr"]

plot(density(alphaNbr_post),
     main="Posterior of Neighbour Effect")
abline(v=0, lty=2)

#Marginal covariate Effects


beta1_post <- samples[,"beta1"]

forest_seq <- seq(0,1,length=100)

psi_curve <- sapply(forest_seq, function(f){
  mean(plogis(mean(samples[,"beta0"]) +
                mean(beta1_post)*f +
                mean(samples[,"beta2"])*0.5))
})

plot(forest_seq, psi_curve, type="l",
     lwd=2,
     xlab="Forest Cover",
     ylab="Occupancy Probability",
     main="Marginal Effect of Forest")

# animated map 

library(gganimate)

df_anim <- data.frame()

for(t in 1:T){
  mat <- to_grid(psi_hat[,t])
  tmp <- melt(mat)
  colnames(tmp) <- c("row","col","value")
  tmp$time <- t
  df_anim <- rbind(df_anim, tmp)
}

p <- ggplot(df_anim, aes(col,row,fill=value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  coord_fixed() +
  scale_y_reverse() +
  transition_states(time) +
  labs(title="Posterior Occupancy - Time {closest_state}")

animate(p)

#posterior predictive check

y_rep_prob <- est_z * plogis(mean(samples[,"theta0"]))

hist(rowMeans(obs_detect),
     main="Observed vs Predicted Detection",
     col="grey")

abline(v=mean(y_rep_prob), col="red", lwd=2)


