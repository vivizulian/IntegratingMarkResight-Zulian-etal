
##############################################################################
####### Code to estimate the Red Snapper abundance on Chicken Rock ###########

## Manuscript:Integrating mark-resight and count data to estimate effective sampling area and fish density.
## Authors: Viviane Zulian, Krishna Pacifici, Nathan M. Bacheler, Jeffrey A. Buckel, William F. Patterson III,
## Brian J. Reich, Kyle W. Shertzer, Nathan J. Hostetter

## Code by: Zulian, V.; Hostetter, N.; Pacifici, K.
## Last edited on July 12, 2023, by V. Zulian.

## Abundance estimates according to two different models: the Standard Binomial N-Mixture model,
## and Marked N-Mixture model

## Last modified on July 5, 2023, by V. Zulian.

## Read the ChickenRockData.RData file ##
# The file includes:
# The number of frames from each period (columns) in each camera (rows) - Nframes
# The site covariates for each camera in each period - chick_rock_covs
# The number of marked individuals detected in each camera, frame and period - m
# The number of unmarked individuals detected in each camera, frame and period - u
# The number of unknown (if marked or unmarked) individuals detected in each camera, frame and period - y
# The total of detections by camera and the number of available individuals around the camera based on VPS data - Ndetec
# The total number of marked individuals available in the entire Chicken Rock by period - Nreef

#load the data set:
#load(".../RedSnapper_ChickenRockData.RData")

library(jagsUI)

#MCP from each individual divided by 1000000 to convert to km2:
#P.S.: We have values for 15 individuals instead of 16, because one of the individuals had only 2 detections. 
MCP <- c(38883,4677,3263,31568,56705,3108,60897,8698,121772,3600,62192,50985,5719,7120,5940)/1000000

## Model 1 - Binomial N-Mixture Model

cat("
  model {

    #Likelihoods
    for (j in 1:J) {		# cameras
      Nlocal[j] ~ dpois(lambda)
      for (k in 1:K) {	# frames
        n[j,k] ~ dbin(p[j], Nlocal[j])
          
      }#K: frames
      logit(p[j]) <- alpha0[curr[j]]
    }#J: cameras

    #Derived abundance for the entire area
    #Not possible to estimate the total abundance, because we don't know ESA of each camera.
    Ntotal ~ dpois(lambda*Area) 

    for(i in 1:15){
      MCP_area[i] ~ dgamma(a, b)
    }

    a ~ dgamma(0.1,0.1) # prior
    b ~ dgamma(0.1,0.1) # prior

    meanMCP <- a/b # expected value
    
    D <- lambda/meanMCP

    # Priors
    for(r in 1:ncurr){
      alpha0[r] <- logit(p_alpha0[r])
      p_alpha0[r] ~ dbeta(1,1)
    }
    lambda ~ dgamma(0.01, 0.01)
    
  }#model", fill=TRUE, file="NMix.txt")

# JAGS data
str(jags.data <- list(J=nrow(n), K=ncol(n), n=n, Area = Area_StateSpace, 
                      curr=curDir+1, ncurr=max(curDir)+1, MCP_area = MCP))

Nst <- apply(n, 1, max, na.rm=T)+3
inits <- function() list(Nlocal = Nst)

# Parameters monitored
params <- c("lambda", "meanMCP", "D", "alpha0", "Nlocal", "Ntotal", "p", "eff.area")

# MCMC settings
#nc <- 3; nt <- 10;  ni <- 500;  nb <- 50;  n.adapt <- 20 #test
nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 

# Call JAGS
library(jagsUI)
out <- jags(jags.data, inits, params, "NMix.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
            n.adapt = n.adapt, parallel=T)

## The sum of Nlocal is not equivalent to the Ntotal, because we don't have information about
## the area sampled of Nlocal.

## Plot the chains:
library(mcmcOutput)
mcmcOutput::diagPlot(out, params = c("lambda", "alpha0", "p")) #Check convergence




## Model 2 - Marked N-Mixture Model

## Model using the marked/unmarked/unknown data sets:
## MODEL Known marked population. Counts record marked, unmarked, and unknown. 
## Estimate proportion of population sampled using known Mtotal

# JAGS data
str(jags.data <- list(J=nrow(n), K=ncol(n), m=mar, u=unm, y=y, 
                      Area_StateSpace = Area_StateSpace, 
                      Mtotal=Mtotal, curDir=curDir+1, ncurr=max(curDir)+1))

cat("
  model {
           
    #Likelihoods
    # processes at state-space level
    Ntotal ~ dpois(lambda*Area_StateSpace) 
    Mtotal ~ dbin(phi, Ntotal)
    Utotal <- Ntotal-Mtotal

    # Abundance at camera sites
    for (j in 1:J) {		
      M[j] ~ dbin(delta[j], Mtotal) 
      U[j] ~ dbin(delta[j], Utotal) 
      N[j] <- M[j]+ U[j]

      # counts by fram
      for (k in 1:K) {	# frames
        m[j,k] ~ dbin(p*theta, M[j])
        u[j,k] ~ dbin(p*theta, U[j])
        y[j,k] ~ dbin(p*(1-theta), N[j])
      }#K occasions
    }#J traps
  
    # Priors
    lambda ~ dgamma(0.01, 0.01) # individuals per unit space
    p ~ dbeta(1,1)              # detection prob given in superpop.
    phi ~ dbeta(1,1)            # prob marked
    theta ~ dbeta(1,1)          # prob marking status is observed if an individual is counted

    # delta
    for(j in 1:J){
      logit(delta[j]) <- gamma0[curDir[j]] + eps[j]
      eps[j] ~ dnorm(0, tau)
    }
       
    tau <- pow(sd, 2)
    sd ~ dgamma(0.01,0.01)
      
    for(r in 1:ncurr){
      gamma0[r] <- logit(p_gamma0[r])
      p_gamma0[r] ~ dbeta(1,1)
    }

  }#model", fill=TRUE, file="MarkedNMix.txt")


#intial values
M_inits <- apply(mar, 1, max, na.rm=T)+2
U_inits <- apply(unm, 1, max, na.rm=T)+2
Nst <- sum(n, na.rm=T)

inits <- function() list(p=runif(1,.06,0.1), 
                         Ntotal=Nst,
                         #lambda=runif(1,50,100),
                         U=U_inits,
                         M=M_inits, 
                         phi=runif(1,.15,0.30), 
                         theta=runif(1,.6,0.9)) #Ntotal=Nst,

# Parameters monitored
params <- c("lambda", "p","theta","phi", "gamma0" ,"Ntotal",
            "delta", "sd", "p_gamma0")#, "M", "U", "N")

#MCMC settings
#nc <- 3; nt <- 5;  ni <- 2000;  nb <- 500;  n.adapt <- 200 test
nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 

# Call JAGS 
out2 <- jags(jags.data, inits, params, "MarkedNMix.txt", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
            n.adapt = n.adapt, parallel=T) 

## Plot the chains:
library(mcmcOutput)
mcmcOutput::diagPlot(out2, params = c("delta", "lambda", "p", "phi", "theta", "gamma0"))


