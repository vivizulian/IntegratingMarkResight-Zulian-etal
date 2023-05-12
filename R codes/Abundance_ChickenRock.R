
#### Code to estimate the Red Snapper abundance on Chicken Rock ####
#### according to two different models: the Standard Binomial N-Mixture model, ####
#### and the modified version of the Multinomial N-Mixture model ####
#### I've used only the sampling period number 5, which is the one with the highest number ####
#### of marked individuals detected. ####
#### For the Binomial N-Mixture model, I've included a Bayesian p-value assessment. ####
#### 
#### Last modified on Dec 7, 2022, by V. Zulian. ####

## Read the ChickenRockData.RData file ##
#The file includes:
#The number of frames from each period (columns) in each camera (rows) - Nframes
#The site covariates for each camera in each period - chick_rock_covs
#The number of marked individuals detected in each camera, frame and period - m
#The number of unmarked individuals detected in each camera, frame and period - u
#The number of unknown (if marked or unmarked) individuals detected in each camera, frame and period - y
#The total of detections by camera and the number of available individuals around the camera based on VPS data - Ndetec
#The total number of marked individuals available in the entire Chicken Rock by period - Nreef

rm(list=ls())

load("Q:/My Drive/Red Snapper/RCodes/Red-Snapper-Abundance/Data/ChickenRockData.RData")

#Run using only the period #5:
#Take subset of the data:
unm <- u[,1:max(Nframes[,5]),5]
mar <- m[,1:max(Nframes[,5]),5]
y <- y[,1:max(Nframes[,5]),5]
#
# #Take a subset of frames and delete camera 13:
unm <- unm[-13,1:24]
mar <- mar[-13,1:24]
y <- y[-13,1:24]
       
n <- matrix(NA, nrow(y), ncol(y))

for(k in 1:nrow(n)){
  for(l in 1:ncol(n)){
    n[k,l] <- y[k,l] + unm[k,l] + mar[k,l]  #total counts
  }
}

Mtotal <- Nreef[5]

#Current direction covariate
curDir <- chic_rock_covs[,"currDirec",5] #period 5
curDir <- curDir[!is.na(curDir)] 

Area_StateSpace <- 1.6
x <- curDir

library(jagsUI)

#Models to test

#####################################
#Model 1 - Binomial N-Mixture Model
#####################################

cat("
    model {

      #Likelihoods
      for (j in 1:J) {		# cameras
        Nlocal[j] ~ dpois(lambda)
        for (k in 1:K) {	# frames
          n[j,k] ~ dbin(p[j], Nlocal[j])
          
        }#K
        logit(p[j]) <- alpha0[curr[j]] #+ eps[j]
        #eps[j] ~ dnorm(0, tau)
      }#J

    #Derived abundance for the entire area
    Ntotal ~ dpois(lambda*Area) 

    # Derived 'effective sampling area'
    #eff.area <- p*Area
    #p.eff.area <- eff.area/Area #estimated proportion of area surveyed by one frame
  
    # Priors
    tau <- pow(sd, 2)
    sd ~ dgamma(0.01,0.01)
      
    for(r in 1:ncurr){
      alpha0[r] ~ dnorm(0, 0.01)
    }
      
    lambda ~ dgamma(0.01, 0.01)
    
  }#model", fill=TRUE, file="Nmix.txt")

# JAGS data
str(jags.data <- list(J=nrow(n), K=ncol(n), n=n, Area = Area_StateSpace, 
                      curr=x+1, ncurr=max(x)+1))

Nst <- apply(n, 1, max, na.rm=T)+3
inits <- function() list(Nlocal = Nst) #lambda=runif(1,.01,0.02)

# Parameters monitored
params <- c("lambda", "alpha0", "sd", "Nlocal", "Ntotal", "p", "eff.area"

#MCMC settings
#nc <- 3; nt <- 10;  ni <- 500;  nb <- 50;  n.adapt <- 20
nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 


# Call JAGS
out <- jags(jags.data, inits, params, "Nmix.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
            n.adapt = n.adapt, parallel=T)

## The sum of Nlocal is not equivalent to the Ntotal, because we don't have information about
## the area sampled of Nlocal.

## Plot the chains:
library(mcmcOutput)
mcmcOutput::diagPlot(out, params = c("lambda", 'alpha1', "alpha0", "p"))

library(here)

#here('Output', "ChickenRock")
save(out, file="ResultNMixt_ChickenRockRandomEff.RData")

############################################
#Model 2 - Multinomial N-Mixture Model
############################################

#B. Model using the marked/unmarked/unknown data sets:
#### MODEL Known marked population. Counts record marked, unmarked, and unknown. 
### Estimate proportion of population sampled using known Mtotal

# JAGS data
str(jags.data <- list(J=nrow(n), K=ncol(n), m=mar, u=unm, y=y, 
                      Area_StateSpace = Area_StateSpace, 
                      Mtotal=Mtotal, x=x+1, ncurr=max(x)+1))

cat("
    model {
           
      #Likelihoods

      # processes at state-space level
      Ntotal ~ dpois(lambda*Area_StateSpace) # is this just a prediction?
      Mtotal ~ dbin(phi, Ntotal) #!!! Mtotal is known!!!
      Utotal <- Ntotal-Mtotal

      # Abundance at camera sites
      for (j in 1:J) {		
        M[j] ~ dbin(delta[j], Mtotal) 
        U[j] ~ dbin(delta[j], Utotal) 
        N[j] <- M[j]+ U[j]

       # counts by frame. lack of independence among count types may affect CrI coverage? 
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

      # options for delta
      # constant
      #delta ~ dbeta(1,1)

      # # parameterization 4: logit-scale effect.   
      for(j in 1:J){
       logit(delta[j]) <- beta0[x[j]] + eps[j]
       eps[j] ~ dnorm(0, tau)
      }
       
      tau <- pow(sd, 2)
      sd ~ dgamma(0.01,0.01)
      
      for(r in 1:ncurr){
        beta0[r] ~ dnorm(0, 0.01)
      }

    }#model", fill=TRUE, file="NmixWithKnownMarked.txt")


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
                         theta=runif(1,.6,0.9)) #beta1= runif(3,-1,0) Ntotal=Nst,

# Parameters monitored
params <- c("lambda", "p","theta","phi", "beta0" ,"Ntotal",
            "delta", "sd", "M", "U", "N")

#MCMC settings
#nc <- 3; nt <- 5;  ni <- 2000;  nb <- 500;  n.adapt <- 200
nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 

# Call JAGS 
out2 <- jags(jags.data, inits, params, "NmixWithKnownMarked.txt", 
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
            n.adapt = n.adapt, parallel=T) 

## Plot the chains:
library(mcmcOutput)
mcmcOutput::diagPlot(out2, params = c("delta", "lambda", "p", "phi", "theta"))

resultsChickenRockPeriod1 <- list(out, out2)
save(resultsChickenRockPeriod1, file="Output/ChickenRock/resultsChickenRockPeriod1.RData")


##################################################################

##_____________CROSS-VALIDATION_______________________________________________##
## Run the model with 80% of the data (training) and predict for the 20% left (testing)
#____________________________________
# Plan is this:
# (1) Randomly exclude 20% of the data - take into NA
# (2) Fit the model and save the posterior samples for the 20% sampled counts
# (3) Repeat n times 
# (4) Plot observed Mean counts vs. predicted Mean counts

#Exclude the counts from one camera each time:
ncam <- 24

#Try to create an empty list to hold the results
resultsMarkedNMixRaEff <- list()

## loop through k number of cameras
for(k in 1:ncam){
  
  # Make a copy of data sets
  marB <- mar
  unmB <- unm
  yB <- y

  # Take one camera out
  marB[k,] <- NA
  unmB[k,] <- NA
  yB[k,] <- NA
  
  # JAGS data
  str(jags.data <- list(J=nrow(n), K=ncol(n), m=marB, u=unmB, y=yB, 
                        Area_StateSpace = Area_StateSpace, 
                        Mtotal=Mtotal, x=x+1, ncurr=max(x)+1))
  
  #intial values
  M_inits <- apply(marB, 1, max, na.rm=T)+2
  M_inits[M_inits == '-Inf'] <- 2
  U_inits <- apply(unmB, 1, max, na.rm=T)+2
  U_inits[U_inits == '-Inf'] <- 2
  Nst <- sum(n, na.rm=T)
  
  inits <- function() list(p=runif(1,.06,0.1), 
                           Ntotal=Nst,
                           #lambda=runif(1,50,100),
                           U=U_inits,
                           M=M_inits, 
                           phi=runif(1,.15,0.30), 
                           theta=runif(1,.6,0.9)) #beta1= runif(3,-1,0) Ntotal=Nst,
  
  # Counter
  cat(paste('\n*** Trial nr.', k, 'out of 24 ***\n'))
  
  # Parameters monitored
  params <- c("lambda", "p","theta","phi", "beta0" ,"Ntotal",
              "delta", "sd", "m", "u", "y")
  
  #MCMC settings
  #nc <- 3; nt <- 5;  ni <- 2000;  nb <- 500;  n.adapt <- 200
  nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 
  
  # Call JAGS 
  fm2 <- jags(jags.data, inits, params, "NmixWithKnownMarked.txt", 
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
               n.adapt = n.adapt, parallel=T) 
  

  # Save predicted Counts for each frame in each MCMC iteration
  #marked individuals
  pred_mar <- apply(fm2$sims.list$m, c(2,3), mean) 
  var_mar <- apply(fm2$sims.list$m, c(2,3), var)
  predMar.CRI025 <- apply(fm2$sims.list$m, c(2,3), function(x) quantile (x, 0.025, na.rm=T))
  predMar.CRI975 <- apply(fm2$sims.list$m, c(2,3), function(x) quantile (x, 0.975, na.rm=T))
  
  #unmarked individuals
  pred_unmar <- apply(fm2$sims.list$u, c(2,3), mean) 
  var_unmar <- apply(fm2$sims.list$u, c(2,3), var)
  predUnmar.CRI025 <- apply(fm2$sims.list$u, c(2,3), function(x) quantile (x, 0.025, na.rm=T))
  predUnmar.CRI975 <- apply(fm2$sims.list$u, c(2,3), function(x) quantile (x, 0.975, na.rm=T))
  
  #unknown individuals
  pred_unk <- apply(fm2$sims.list$y, c(2,3), mean) 
  var_unk <- apply(fm2$sims.list$y, c(2,3), var)
  predUnk.CRI025 <- apply(fm2$sims.list$y, c(2,3), function(x) quantile (x, 0.025, na.rm=T))
  predUnk.CRI975 <- apply(fm2$sims.list$y, c(2,3), function(x) quantile (x, 0.975, na.rm=T))
  
  tmp <- list(pred_mar=pred_mar, var_mar=var_mar, predMar.CRI025=predMar.CRI025, predMar.CRI975=predMar.CRI975,
              pred_unmar=pred_unmar, var_unmar=var_unmar, predUnmar.CRI025=predUnmar.CRI025, predUnmar.CRI975=predUnmar.CRI975,
              pred_unk=pred_unk, var_unk=var_unk, predUnk.CRI025=predUnk.CRI025, predUnk.CRI975=predUnk.CRI975,
              fm2 = fm2)
  
  name <- paste('trial:', k)
  resultsMarkedNMixRaEff[[name]] <- tmp
}





#away=0, sideways=1, and towards=2

away <- as.data.frame(cbind(out2$sims.list$beta0[,1], rep("away", length(out2$sims.list$beta0[,1]))))
sideways <- as.data.frame(cbind(out2$sims.list$beta0[,2], rep("sideways", length(out2$sims.list$beta0[,2]))))
towards <- as.data.frame(cbind(out2$sims.list$beta0[,3], rep("towards", length(out2$sims.list$beta0[,3]))))

current <- rbind(away, sideways, towards)

current$V3 <- plogis(as.numeric(current$V1))
current$area <- current$V3*Area_StateSpace


#Plot the results:
require(ggplot2)
library(ggforce)
(Curr <- ggplot(current, aes(y = as.numeric(area), x = V2, fill=V2)) + 
    geom_violin() + theme_bw() +
    #geom_sina(alpha=0.01, size=0.01) + 
    geom_boxplot(width=0.04, fill="white", outlier.shape = NA) + 
    #stat_summary(fun = "median", geom = "crossbar", width = 0.1, colour = "black", position = position_dodge(0.9)) +
    #geom_line(aes(y = meanNstarSigma2, x=covs), size=1) + 
    #geom_hline(yintercept=c(Nstar_Sigma0, Nstar_Sigma1, meanNstarSigma2), linetype="dashed", size=1) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), text = element_text(size = 17), 
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          legend.position = NULL) + #c(0.15, 0.87)) + 
    ylim(0,1) +
    labs(y = "Effective sampling area (km2)", x="Current direction") + 
    scale_fill_manual(values = c("#CCCCCC", "#CCCCCC", "#CCCCCC")) +  
    scale_x_discrete(labels = c("Away", "Sideways", "Towards")))


PerCam <- as.data.frame(cbind(out2$mean$delta, out2$mean$delta*Area_StateSpace, curDir))


#Plot the results:
require(ggplot2)
library(ggforce)
(Curr <- ggplot(PerCam, aes(y = as.numeric(V2), x = as.factor(curDir), fill=as.factor(curDir))) + 
    geom_violin() + theme_bw() +
    geom_sina(alpha=1, size=0.9) + 
    geom_boxplot(width=0.05, fill="white", outlier.shape = NA) + 
    #stat_summary(fun = "median", geom = "crossbar", width = 0.1, colour = "black", position = position_dodge(0.9)) +
    #geom_line(aes(y = meanNstarSigma2, x=covs), size=1) + 
    #geom_hline(yintercept=c(Nstar_Sigma0, Nstar_Sigma1, meanNstarSigma2), linetype="dashed", size=1) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), text = element_text(size = 17), 
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          legend.position = NULL) + #c(0.15, 0.87)) + 
    ylim(0,1.6) +
    labs(y = "Effective sampling area (km2)", x="Current direction") + 
    scale_fill_manual(values = c("#CCCCCC", "#CCCCCC", "#CCCCCC")) +  
    scale_x_discrete(labels = c("Away", "Sideways", "Towards")))


