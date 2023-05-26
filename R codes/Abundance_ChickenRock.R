
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

#load the data set:
#load(".../RedSnapper_ChickenRockData.RData")

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
        logit(p[j]) <- alpha0[curr[j]]
      }#J

    #Derived abundance for the entire area
    Ntotal ~ dpois(lambda*Area) 

    # Derived 'effective sampling area'
    #eff.area <- p*Area
    #p.eff.area <- eff.area/Area #estimated proportion of area surveyed by one frame
  
    # Priors
    for(r in 1:ncurr){
      alpha0[r] ~ dnorm(0, 0.01)
    }
      
    lambda ~ dgamma(0.01, 0.01)
    
  }#model", fill=TRUE, file="NMix.txt")

# JAGS data
str(jags.data <- list(J=nrow(n), K=ncol(n), n=n, Area = Area_StateSpace, 
                      curr=curDir+1, ncurr=max(curDir)+1))

Nst <- apply(n, 1, max, na.rm=T)+3
inits <- function() list(Nlocal = Nst) #lambda=runif(1,.01,0.02)

# Parameters monitored
params <- c("lambda", "alpha0", "Nlocal", "Ntotal", "p", "eff.area")

#MCMC settings
#nc <- 3; nt <- 10;  ni <- 500;  nb <- 50;  n.adapt <- 20 #test
nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 


# Call JAGS
out <- jags(jags.data, inits, params, "NMix.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
            n.adapt = n.adapt, parallel=T)

## The sum of Nlocal is not equivalent to the Ntotal, because we don't have information about
## the area sampled of Nlocal.

## Plot the chains:
library(mcmcOutput)
mcmcOutput::diagPlot(out, params = c("lambda", "alpha0", "p"))


############################################
#Model 2 - Marked N-Mixture Model
############################################

#B. Model using the marked/unmarked/unknown data sets:
#### MODEL Known marked population. Counts record marked, unmarked, and unknown. 
### Estimate proportion of population sampled using known Mtotal

# JAGS data
str(jags.data <- list(J=nrow(n), K=ncol(n), m=mar, u=unm, y=y, 
                      Area_StateSpace = Area_StateSpace, 
                      Mtotal=Mtotal, curDir=curDir+1, ncurr=max(curDir)+1))

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
ncam <- 31

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
  params <- c("lambda", "p","theta","phi", "gamma0" ,"Ntotal",
              "delta", "sd", "p_gamma0") #, "m", "u", "y")
  
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


#### Cross-Validation processing ####

### Random effect model:

#load results
load("Q:/My Drive/Red Snapper/RCodes/Red-Snapper-Abundance/Output/ChickenRock/resultsMarkedNMixRaEff_crossVal.RData")

pred_mar <- var_mar <- predMar.CRI025 <- predMar.CRI975 <- pred_unmar <- var_unmar <- predUnmar.CRI025 <- predUnmar.CRI975 <-
  pred_unk <- var_unk <- predUnk.CRI025 <- predUnk.CRI975 <- matrix(NA, 31, 24)

ncam = 24
for(i in 1:ncam){
  pred_mar[i,] <- resultsMarkedNMixRaEff[[i]]$pred_mar[i,]
  var_mar[i,] <- resultsMarkedNMixRaEff[[i]]$var_mar[i,]
  predMar.CRI025[i,] <- resultsMarkedNMixRaEff[[i]]$predMar.CRI025[i,]
  predMar.CRI975[i,] <- resultsMarkedNMixRaEff[[i]]$predMar.CRI975[i,]
  pred_unmar[i,] <- resultsMarkedNMixRaEff[[i]]$pred_unmar[i,]
  var_unmar[i,] <- resultsMarkedNMixRaEff[[i]]$var_unmar[i,]
  predUnmar.CRI025[i,] <- resultsMarkedNMixRaEff[[i]]$predUnmar.CRI025[i,]
  predUnmar.CRI975[i,] <- resultsMarkedNMixRaEff[[i]]$predUnmar.CRI975[i,]
  pred_unk[i,] <- resultsMarkedNMixRaEff[[i]]$pred_unk[i,]
  var_unk[i,] <- resultsMarkedNMixRaEff[[i]]$var_unk[i,]
  predUnk.CRI025[i,] <- resultsMarkedNMixRaEff[[i]]$predUnk.CRI025[i,]
  predUnk.CRI975[i,] <- resultsMarkedNMixRaEff[[i]]$predUnk.CRI975[i,]
}



















mean(out2$sims.list$lambda)
mean(out2$sims.list$Ntotal)

median(out2$sims.list$p)
quantile(out2$sims.list$p, probs=c(0.025, 0.975))

median(out2$sims.list$theta)
quantile(out2$sims.list$theta, probs=c(0.025, 0.975))

median(out2$sims.list$phi)
quantile(out2$sims.list$phi, probs=c(0.025, 0.975))

#away=0, sideways=1, and towards=2

away <- as.data.frame(cbind(out2$sims.list$p_gamma0[,1], rep("away", length(out2$sims.list$p_gamma0[,1]))))
round(mean(as.numeric(away$V1)*Area_StateSpace),3)
round(quantile(as.numeric(away$V1)*Area_StateSpace, probs=c(0.025, 0.975)),3)

sideways <- as.data.frame(cbind(out2$sims.list$p_gamma0[,2], rep("sideways", length(out2$sims.list$p_gamma0[,2]))))
round(mean(as.numeric(sideways$V1)*Area_StateSpace),3)
round(quantile(as.numeric(sideways$V1)*Area_StateSpace, probs=c(0.025, 0.975)),3)

towards <- as.data.frame(cbind(out2$sims.list$p_gamma0[,3], rep("towards", length(out2$sims.list$p_gamma0[,3]))))
round(mean(as.numeric(towards$V1)*Area_StateSpace),3)
round(quantile(as.numeric(towards$V1)*Area_StateSpace, probs=c(0.025, 0.975)),3)

current <- rbind(away, sideways, towards)

Area_StateSpace <- 1.6
#current$V3 <- plogis(as.numeric(current$V1))
current$area <- as.numeric(current$V1)*Area_StateSpace

library(tidyverse)
ci <- current %>%
  group_by(V2) %>%
  summarise(lower.ci = round(quantile(area, probs = c(0.025)),3),
            upper.ci = round(quantile(area, probs = c(0.975)),3),
            med = round(median(area),3))


#Plot the results:
require(ggplot2)
library(ggforce)
(Curr <- ggplot() + 
    geom_violin(data = current, aes(y = as.numeric(area), x = V2, fill=V2), color=NA) + theme_bw() +
    #stat_summary(fun = "median", geom = "pointrange", position = position_dodge(0.9)) + # fun.min = min, fun.max = max, linewidth = 2, size = 0.5, colour = "black") + #geom = "pointrange", , position = position_dodge(0.9)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), text = element_text(size = 17), 
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          legend.position = NULL) + #c(0.15, 0.87)) + 
    ylim(0,0.4) +
    geom_linerange(data = ci, aes(x = V2, ymax = upper.ci, ymin = lower.ci), color="black", lwd = 1) +
  geom_point(data = ci, aes(x = V2, y = med), color = 'black', size = 3) +
    labs(y = "Effective sampling area (km2)", x="Current direction") + 
    scale_fill_manual(values = c("#CCCCCC", "#CCCCCC", "#CCCCCC")) +  
    scale_x_discrete(labels = c("Away", "Sideways", "Towards")))

#
#geom_sina(alpha=0.01, size=0.01) + 
#geom_boxplot(width=0.04, fill="white", outlier.shape = NA) + 
#geom_point(current, mapping = aes(as.numeric(area), )) +

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




## sims.list of delta is all the posteriors for the 31 cameras

ESAcameras <- out2$sims.list$delta*Area_StateSpace
meanESAcam <- colMeans(ESAcameras)
medianESAcam <- apply(ESAcameras, 2, median)
Ci25ESAcam <- apply(ESAcameras, 2, quantile, probs=c(0.025))
Ci975ESAcam <- apply(ESAcameras, 2, quantile, probs=c(0.975))

ESAcams <- cbind(seq(1:31), meanESAcam, medianESAcam, Ci25ESAcam, Ci975ESAcam, curDir)
ESAcams[,6] <- ifelse(ESAcams[,6]==0, paste("Away"), ESAcams[,6])
ESAcams[,6] <- ifelse(ESAcams[,6]==1, paste("Sideways"), ESAcams[,6])
ESAcams[,6] <- ifelse(ESAcams[,6]==2, paste("Towards"), ESAcams[,6])

ESAcams <- as.data.frame(ESAcams)

library(tidyverse)
ci <- ESAcams %>%
  group_by(curDir) %>%
  summarise(lower.ci = round(quantile(as.numeric(meanESAcam), probs = c(0.025)),3),
            upper.ci = round(quantile(as.numeric(meanESAcam), probs = c(0.975)),3),
            med = round(mean(as.numeric(meanESAcam)),3))

ci$med

#Plot the results:
require(ggplot2)
library(ggforce)
(Curr <- ggplot(ESAcams, aes(y = as.numeric(meanESAcam), x = as.factor(V1), color=curDir)) + 
    
    #annotate(geom = "rect", ymin=ci$lower.ci[1], ymax=ci$upper.ci[1],
    #         xmin=-Inf, xmax=Inf, fill="gray", alpha=0.4) +
    #annotate(geom = "rect", ymin=ci$lower.ci[2], ymax=ci$upper.ci[2],
    #         xmin=-Inf, xmax=Inf, fill="orange", alpha=0.1) +
    #annotate(geom = "rect", ymin=ci$lower.ci[3], ymax=ci$upper.ci[3],
    #         xmin=-Inf, xmax=Inf, fill="blue", alpha=0.1) +
    
    geom_point(size=3) + theme_bw() +
    geom_linerange(data = ESAcams, aes(x = V1, ymax = as.numeric(Ci975ESAcam), ymin = as.numeric(Ci25ESAcam)),lwd = 1) +
    geom_hline(yintercept=c(ci$med[1], ci$med[2], ci$med[3]), linetype="dashed", size=1, color=c("#000000", "#E69F00", "#56B4E9")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          plot.title = element_text(hjust = 0.5), text = element_text(size = 17), 
          axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
          legend.position = c(0.15, 0.87)) + 
    ylim(0,1.6) +
    labs(y = "Effective sampling area (km2)", x="Camera", colour="Current direction") + 
    scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")))# +  
    #scale_x_discrete(labels = c("Away", "Sideways", "Towards")))
