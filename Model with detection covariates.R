
  
# ==========================================================

# FULL DYNAMIC OCCUPANCY + MULTI-DETECTION MODEL

# ==========================================================


############################################################
## 1. SETUP
############################################################

library(nimble)
library(coda)
library(ggplot2)

set.seed(123)

nrow_grid <- 20      # reduced for speed
ncol_grid <- 20
N <- nrow_grid * ncol_grid
T <- 4

nScout <- 3
nDrone <- 3
nWeek  <- 8
maxCam <- 2


  # 2. BUILD ADJACENCY MATRIX
  
build_adjacency <- function(nrow, ncol){
  N <- nrow * ncol
  W <- matrix(0, N, N)
  index <- function(r,c) (r-1)*ncol + c
  
  for(r in 1:nrow){
    for(c in 1:ncol){
      i <- index(r,c)
      neighbors <- list(
        c(r-1,c),c(r+1,c),c(r,c-1),c(r,c+1),
        c(r-1,c-1),c(r+1,c-1),c(r+1,c+1),c(r-1,c+1)
      )
      for(nb in neighbors){
        rr <- nb[1]; cc <- nb[2]
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


  # 3. SIMULATE ECOLOGICAL COVARIATES
  
forest    <- runif(N)
distCrop  <- runif(N)
control   <- matrix(rbinom(N*T,1,0.8), N, T)

  # 4. SIMULATE DETECTION COVARIATES
  
scoutCond <- array(sample(1:3, N*T*nScout, replace=TRUE),
                   c(N,T,nScout))

droneCond <- array(sample(1:3, N*T*nDrone, replace=TRUE),
                   c(N,T,nDrone))

camFOV <- array(sample(1:3, N*T*maxCam*nWeek, replace=TRUE),
                c(N,T,maxCam,nWeek))

  # 5. TRUE PARAMETERS
  
true <- list(
  beta0=-1, beta1=1.2, beta2=0.3,
  alpha0=-2, alphaNbr=0.6, alpha2=0.5, alpha3=-0.3,
  delta0=-1, delta1=1.2,
  
  theta0=-1, theta_adequate=0.4, theta_good=0.8,
  eta0=-1, eta_adequate=0.5, eta_good=1,
  phi0=-1, phi_old=0.5, phi_fresh=1,
  sigma_cam=0.4
)

  # 6. SIMULATE LATENT OCCUPANCY
  
z <- matrix(0,N,T)

psi1 <- plogis(true$beta0 +
                 true$beta1*forest +
                 true$beta2*distCrop)

z[,1] <- rbinom(N,1,psi1)

for(t in 2:T){
  nOccNbr <- W %*% z[,t-1]
  
  gamma <- plogis(true$alpha0 +
                    true$alphaNbr*nOccNbr +
                    true$alpha2*forest +
                    true$alpha3*distCrop)
  
  epsilon <- plogis(true$delta0 +
                      true$delta1*control[,t])
  
  psi <- z[,t-1]*(1-epsilon) +
    (1-z[,t-1])*gamma
  
  z[,t] <- rbinom(N,1,psi)
}

  # 7. SIMULATE DETECTION
  
  ### Scout
  
yScout <- array(0,c(N,T,nScout))
for(i in 1:N){
  for(t in 1:T){
    for(j in 1:nScout){
      cond <- scoutCond[i,t,j]
      lp <- true$theta0 +
        (cond==2)*true$theta_adequate +
        (cond==3)*true$theta_good
      yScout[i,t,j] <- rbinom(1,1,z[i,t]*plogis(lp))
    }
  }
}

  ### Drone
  
yDrone <- array(0,c(N,T,nDrone))
for(i in 1:N){
  for(t in 1:T){
    for(j in 1:nDrone){
      cond <- droneCond[i,t,j]
      lp <- true$eta0 +
        (cond==2)*true$eta_adequate +
        (cond==3)*true$eta_good
      yDrone[i,t,j] <- rbinom(1,1,z[i,t]*plogis(lp))
    }
  }
}

  ### Camera (with random effect)
  
camRE_true <- matrix(rnorm(N*maxCam,0,true$sigma_cam),
                     N,maxCam)

yCam <- array(0,c(N,T,maxCam,nWeek))

for(i in 1:N){
  for(t in 1:T){
    for(k in 1:maxCam){
      for(w in 1:nWeek){
        fov <- camFOV[i,t,k,w]
        lp <- true$phi0 +
          (fov==2)*true$phi_old +
          (fov==3)*true$phi_fresh +
          camRE_true[i,k]
        yCam[i,t,k,w] <- rbinom(1,1,z[i,t]*plogis(lp))
      }
    }
  }
}

  # 8. NIMBLE MODEL
  
code <- nimbleCode({
  
  # PRIORS
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  
  alpha0 ~ dnorm(0,0.01)
  alphaNbr ~ dnorm(0,0.01)
  alpha2 ~ dnorm(0,0.01)
  alpha3 ~ dnorm(0,0.01)
  
  delta0 ~ dnorm(0,0.01)
  delta1 ~ dnorm(0,0.01)
  
  theta0 ~ dnorm(0,0.01)
  theta_adequate ~ dnorm(0,0.01)
  theta_good ~ dnorm(0,0.01)
  
  eta0 ~ dnorm(0,0.01)
  eta_adequate ~ dnorm(0,0.01)
  eta_good ~ dnorm(0,0.01)
  
  phi0 ~ dnorm(0,0.01)
  phi_old ~ dnorm(0,0.01)
  phi_fresh ~ dnorm(0,0.01)
  
  sigma_cam ~ dunif(0,5)
  tau_cam <- pow(sigma_cam,-2)
  
  for(i in 1:N){
    
    logit(psi[i,1]) <- beta0 +
      beta1*forest[i] +
      beta2*distCrop[i]
    
    z[i,1] ~ dbern(psi[i,1])
    
    for(t in 2:T){
      
      nOccNbr[i,t-1] <- inprod(W[i,1:N], z[1:N,t-1])
      
      logit(gamma[i,t]) <- alpha0 +
        alphaNbr*nOccNbr[i,t-1] +
        alpha2*forest[i] +
        alpha3*distCrop[i]
      
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
      
      for(j in 1:nDrone){
        logit(pD[i,t,j]) <- eta0 +
          equals(droneCond[i,t,j],2)*eta_adequate +
          equals(droneCond[i,t,j],3)*eta_good
        yDrone[i,t,j] ~ dbern(z[i,t]*pD[i,t,j])
      }
      
      for(k in 1:maxCam){
        
        camRE[i,k] ~ dnorm(0,tau_cam)
        
        for(w in 1:nWeek){
          logit(pC[i,t,k,w]) <- phi0 +
            equals(camFOV[i,t,k,w],2)*phi_old +
            equals(camFOV[i,t,k,w],3)*phi_fresh +
            camRE[i,k]
          
          yCam[i,t,k,w] ~ dbern(z[i,t]*pC[i,t,k,w])
        }
      }
    }
  }
})


  # 9. BUILD + RUN MCMC
  
constants <- list(
  N=N,T=T,
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

inits <- list(z=matrix(1,N,T))

model <- nimbleModel(code,constants,data,inits)
cmodel <- compileNimble(model)

conf <- configureMCMC(model,
                      monitors=c("beta0","beta1","beta2",
                                 "alphaNbr","delta1",
                                 "theta0","eta0","phi0",
                                 "sigma_cam"))

mcmc <- buildMCMC(conf)
cmcmc <- compileNimble(mcmc,project=model)

samples <- runMCMC(cmcmc,
                   niter=4000,
                   nburnin=1500,
                   nchains=2,
                   samplesAsCodaMCMC=TRUE)

  # 10. CHECK MIXING
  
plot(samples)
gelman.diag(samples)
effectiveSize(samples)

  # 11. POSTERIOR SUMMARIES
  
summary(samples)

post_means <- summary(samples)$statistics[,1]
print(post_means)


  # 12. VISUALIZE TRUE vs ESTIMATED
  
truth_vec <- c(
  true$beta0,
  true$beta1,
  true$beta2,
  true$alphaNbr,
  true$delta1,
  true$theta0,
  true$eta0,
  true$phi0,
  true$sigma_cam
)

names(truth_vec) <- names(post_means)

df <- data.frame(
  Parameter=names(post_means),
  Posterior=post_means,
  Truth=truth_vec
)

ggplot(df,aes(x=Truth,y=Posterior))+
  geom_point(size=3)+
  geom_abline(slope=1,intercept=0,col="red")+
  theme_bw()+
  ggtitle("Posterior Mean vs True Value")

  # 13. OCCUPANCY MAP VISUALIZATION
  
z_est <- model$z

par(mfrow=c(1,2))
image(matrix(z[,T],nrow_grid,ncol_grid),
      main="True occupancy final")
image(matrix(z_est[,T],nrow_grid,ncol_grid),
      main="Estimated occupancy final")


