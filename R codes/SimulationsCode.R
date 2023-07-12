
##############################################################################
############## Code for simulating and analyzing count data ##################

## Manuscript:Integrating mark-resight and count data to estimate effective sampling area and fish density.
## Authors: Viviane Zulian, Krishna Pacifici, Nathan M. Bacheler, Jeffrey A. Buckel, William F. Patterson III,
## Brian J. Reich, Kyle W. Shertzer, Nathan J. Hostetter

## Code by: Zulian, V.; Hostetter, N.; Pacifici, K.
## Last edited on July 12, 2023, by V. Zulian.

## This code simulates a 100 x 100 state space, with individuals randomly
## distributed under a density (λ) of 0.15 individuals per unit area,
## and a probability of being marked (φ) of 0.1.
## The sampling process consist of 49 sampling sites (J), with four occasions (K) per site.
## Each site had a fixed survey radius (ESA) of 50 units of area. 

## The divided in two simulation approaches: 
## APPROACH 1: data sets simulated under the assumption of no movement 
## of individuals between occasions. 

## APPROACH 2: data sets simulated under the assumption of movement 
## of individuals between occasions: in each occasions, individual 
## movements followed a Bivariate Normal distribution, with standard deviation σ (sigma)

## Simulated data from both approaches is analyzed using a standard N-Mixture model and
## the Marked N-Mixture model.

##________________________________________________________________________________##
## APPROACH 1: data sets simulated under the assumption of no movement 
## of individuals between occasions. 

## Some notes
# S = state space
# J = sites (cameras)
# K = occasions (frames)
# Ntotal = abundance in the state space
# N[j] = abundance at site j
# ylim = xlim = boundaries of S
# Area_StateSpace = area of state space
# area_camera[k] = surveyed area of camera k

## Set a seed for the simulation
set.seed(42+sim)
                          
# Data simulation
lambda <- 0.15	  # individuals per unit space 
phi <- 0.1        # prob of that an individual is marked
p <- 0.5          # 0.5 for models A-C: detection prob | in camera field of view. 
#p <- 0.25        # 0.25 for models D-F: detection prob | in camera field of view. 
theta <- 1        # 1 for model A: prob that marking status is observed
#theta <- 0.75    # 0.75 for models B-F: prob that marking status is observed
                          
# Sampling via cameras
J <- 49 	 # number of sites (cameras)
K <- 4		 # number of occasions (frames)
                          
# State-space
xlim <- c(0,100)	# x-coordinates for reef
ylim <- c(0,100)	# y-coordinates for reef
Area_StateSpace <- diff(xlim)* diff(ylim)
                          
# Camera locations: keep it sufficiently buffered relative to state-space
xy_cam <- expand.grid(seq(20, 80, length=sqrt(J)), seq(20, 80, length=sqrt(J)))
                          
# Area surveyed - known now, set beta1 to 0 and give delta as data to the model.
x <- rnorm(J, 0, 1)              # standardized camera covariate (e.g., turbidity)
beta0 <- qlogis(0.005)           # proportion of state-space sampled by a single camera when x=0
beta1 <- 0 #-0.50                # covariate relationship - if needed
delta <- plogis(beta0 + beta1*x) # delta varies by x if beta1 does not equal 0
area_camera <- delta*Area_StateSpace # area sampled by camera k is now derived
radius <- sqrt(area_camera/pi)
                       
### Generate stochastic realization
n <- m <- u <- y <- matrix(NA, J, K)
                          
# STATE-SPACE: Abundance and marking process
Ntotal <- rpois(1,lambda*Area_StateSpace) # total abundance on the reef
Mtotal <- rbinom(1, Ntotal, phi)          # number of marked individuals on the reef
Utotal <- Ntotal - Mtotal                 # number of unmarked individuals
                          
# CAMERA-LEVEL abundance 
M <- rbinom(J, Mtotal, delta) 
U <- rbinom(J, Utotal, delta) 
N <- M + U
                         
# OBSERVATIONS by occasions (frames)
for(j in 1:J){   # camera
  for(k in 1:K){  # occasion (frame)
    # marked individuals observed
    m_tmp <- rbinom(1, M[j], p) 
    # recorded assignments of observed marked individuals
    m[j,k] <- rbinom(1, m_tmp, theta)  # observed and recorded as marked
    m_unk <- m_tmp - m[j,k]            # observed but marking status was unknown
                              
    # unmarked individuals observed
    u_tmp <- rbinom(1, U[j], p) 
    # recorded assignments of observed unmarked individuals
    u[j,k] <- rbinom(1, u_tmp, theta)  # observed and recorded as unmarked
    u_unk <- u_tmp - u[j,k]            # observed but marking status was unknown
                              
    # unknown observations are derived and a mix of marked and unmarked individuals
    y[j,k] <- m_unk + u_unk            # observed unknown
                              
    # total number of individuals observed
    n[j,k] <- m[j,k] + u[j,k] + y[j,k] 
  } #k
} #j
                          
mean_n <- apply(n, 1, mean)
simsData <- list(N, mean_n) # save the simulated data


### MODEL 1: N-MIXTURE model: ###
cat("
   model {
   
    #Likelihoods
    for (j in 1:J) {		# cameras
      Nlocal[j] ~ dpois(lambda*area[j])
      for (k in 1:K) {	# frames
        n[j,k] ~ dbin(p, Nlocal[j])
      }#K
    }#J
                    
    #Derived abundance for the entire area
    Ntotal ~ dpois(lambda*Area)
                    
    # Priors
    lambda ~ dgamma(0.01, 0.01)
    p ~ dbeta(1,1)
                    
  }#model", fill=TRUE, file="Nmix.txt")
                          
# JAGS data
jags.data <- list(K=K, J=J, n=n, Area=Area_StateSpace, area=area_camera)
Nst <- sum(n)
inits <- function() list(p=runif(1,.6,0.8), Ntotal=Nst)
                          
# Parameters monitored
params <- c("p", "lambda", "Nlocal")
                          
# MCMC settings
#nc <- 3; nt <- 5;  ni <- 200;  nb <- 50;  n.adapt <- 20 # test
nc <- 3; nt <- 50;  ni <- 200000;  nb <- 50000;  n.adapt <- 20000
                          
# Call JAGS
library(jagsUI)
out <- jags(jags.data, inits, params, "Nmix.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
            n.adapt = n.adapt, parallel=T)
                          


### MODEL 2: MARKED N-MIXTURE model: ###   
cat("
  model {
                                 
    # Likelihoods
    # processes at state-space level
    Ntotal ~ dpois(lambda*Area_StateSpace)
    Mtotal ~ dbin(phi, Ntotal) # Mtotal is known
    Utotal <- Ntotal-Mtotal
                      
    # Abundance at camera sites
    for (j in 1:J) {		
      M[j] ~ dbin(delta[j], Mtotal) 
      U[j] ~ dbin(delta[j], Utotal) 
      N[j] <- M[j]+ U[j]
                      
      # Counts by frame k
      for (k in 1:K) {	# frames
        m[j,k] ~ dbin(p*theta, M[j])
        u[j,k] ~ dbin(p*theta, U[j])
        y[j,k] ~ dbin(p*(1-theta), N[j])
      } # K occasions or frames
    } # J sites or cameras 
                            
    # Priors
    lambda ~ dgamma(0.01, 0.01)  # individuals per unit space
    p ~ dbeta(1,1)               # detection prob given in superpop.
    phi ~ dbeta(1,1)             # prob marked
    # theta ~ dbeta(1,1)         # prob marking status is observed if an individual is counted - estimated on models B-F
                      
    # options for delta
    # constant - if no covariates
    # delta ~ dbeta(1,1)
                      
    # parameterization in logit-scale random effect - if covariates
    # for(j in 1:J){
      # logit(delta[j]) <- beta0 + beta1*x[j]       # proportion of population sampled
    # }
    # beta0 <- logit(mean.delta)
    # mean.delta ~ dbeta(1,1)
    # beta1 ~ dnorm(0,.1)
                      
    ## Derive some metrics
    # mean effective area sampled by a camera. if geographically open, effectiveArea will be > area 
    effectiveArea <- delta*Area_StateSpace
                      
    # superpopulation exposed to sampling at a camera
    realized_N_star <- mean(N[1:J])
    expected_N_star <- lambda*Area_StateSpace*delta
                      
  } #model", fill=TRUE, file="NMarked.txt")
                                              
# JAGS data
jags.data <- list(K=K, J=J, m=m, u=u, y=y, Area_StateSpace = Area_StateSpace, 
                  Mtotal=Mtotal, x=x, delta=delta, theta=1)
                          
# Intial values
M_inits <- apply(m, 1, max)+2
U_inits <- apply(u, 1, max)+2
                          
inits <- function() list(p=runif(1,.6,0.8), 
                         lambda=runif(1,.10,0.30),
                         U=U_inits,
                         M=M_inits,
                         phi=runif(1,.15,0.30))#, 
                         # theta=runif(1,.6,0.95))
                         # mean.delta=runif(1,.0001,.001),    # proportion sampled by a camera is very low
                         # beta1= runif(1,-1,0)
 
# MCMC settings                    
# nc <- 3; nt <- 5;  ni <- 200;  nb <- 50;  n.adapt <- 20 # test     
nc <- 3; n.adapt <- 70000; nb <- 5000; ni <- 500000+nb; nt <- 100 
                          
# Parameters monitored
params <- c("lambda", "p","phi", #"mean.delta", "beta1",
            "effectiveArea", "realized_N_star", "expected_N_star" )
                          
# Call JAGS 
out2 <- jags(jags.data, inits, params, "NMarked.txt", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
             n.adapt = n.adapt, parallel=T) 
                          


##________________________________________________________________________________##
## APPROACH 2: data sets simulated under the assumption of movement of individuals
## between occasions. 

## Some notes
# S = statespace
# J = sites (cameras) 
# K = occasions (frames)
# Ntotal = abundance in state space
# N[j] = abundance at site j
# i = 1, 2, ... Ntotal individuals
# s[i,1:2] = activity center of each individual
# u[i,1:2,k] = location during survey k
# ylim = xlim = boundaries of S
# Area_StateSpace = area of state space
# area_camera[k] = surveyed area of camera k

## Set a seed for the simulation
set.seed(42+sim)

# Data simulation 
lambda <- 0.15	  # individuals per unit space 
phi <- 0.1        # prob of that an individual is marked
p <- 0.50         # detection prob | in camera field of view. 
theta <- 0.75     # prob that marking status is observed
sigma <- 2        # 0, 1, or 2 (when sigma = 0, model output is the same as model B from approach 1)

# Sampling via cameras
J <- 49 	# number of sites (cameras)
K <- 4		# number of occasions (frames)

# State-space
xlim <- c(0,100)	# x-coordinates for state space (reef)
ylim <- c(0,100)	# y-coordinates for state space (reef)
Area_StateSpace <- diff(xlim)* diff(ylim)

# Camera locations - keep sufficiently buffered relative to state-space
xy_cam <- expand.grid(seq(20, 80, length=sqrt(J)), seq(20, 80, length=sqrt(J)))

area_camera<- rep(50,J) #area sampled by each camera.
delta <- area_camera/Area_StateSpace # proportion of area sampled by camera k
radius <- sqrt(area_camera/pi)

### Generate stochastic realization
# STATE-SPACE: Abundance and marking process
Ntotal <- rpois(1,lambda*Area_StateSpace) # total abundance on the reef
Mtotal <- rbinom(1, Ntotal, phi)          # number of marked individuals on the reef
Utotal <- Ntotal - Mtotal                 # number of unmarked individuals on the reef

# INDIVIDUAL-LEVEL locations, space use, and interactions with cameras 
# Randomly mark Mtotal individuals
marked <- sample(1:Ntotal, Mtotal) 
mark <- rep(0,Ntotal)
mark[marked]=1      # assign which are the marked fish

# Data placeholders for individual level encounter histories
u_i <- m_i <- y_i <- i_local <- array(NA, c(Ntotal, J, K))

# Activity centers of each individual
s <- cbind(runif(Ntotal, xlim[1], xlim[2]), runif(Ntotal, ylim[1], ylim[2])) 

# Location of each individual during a survey placeholder
mu <- array(NA, c(Ntotal, 2, K)) 

# Simulate realizations
for(i in 1:Ntotal){ # individuals
  mu[i, 1:2, 1:K] <- rnorm(K*2, s[i, 1:2], sd=sigma) # location during a survey
  
  for(k in 1:K){ # occasions
    d <- sqrt(((mu[i,1,k]-xy_cam[1:J,1])^2) + ((mu[i,2,k]-xy_cam[1:J,2])^2)) # distance to traps
    
    # is individual i within the sampling radius of cameras 1:J on survey k
    i_local[i, 1:J, k] <- (d[1:J]<radius[1:J])*1 
    
    # detection can only occur if individual is in field of view
    p_real <- p*i_local[i, 1:J, k]
    
    # detection events
    detected <- rbinom(J,1,p_real)
    
    for(j in 1:J){ #cameras
      u_i[i,j,k] <- rbinom(1,detected[j],theta)*(1-mark[i])  #observed as unmarked | unmarked and observed
      m_i[i,j,k] <- rbinom(1,detected[j],theta)*(mark[i])    #observed as marked | marked and observed
      y_i[i,j,k] <- detected[j] - u_i[i,j,k] - m_i[i,j,k]    #observed, but marking status unknown
    }#J
  }#K
}#i

N <- apply(i_local, c(2,3), sum)                     # number of individuals present during each survey
N_star <- apply(apply(i_local, c(1,2), max), 2, sum) # cumulative number of individuals that used a site

## Observed number of individuals (counts) by camera and occasion 
u <- apply(u_i, c(2,3), sum)     # observed as unmarked
m <- apply(m_i, c(2,3), sum)     # observed as marked
y <- apply(y_i, c(2,3), sum)     # observed as unknown

# Total count (for standard N-mix)
n <- matrix(NA, J, K)
for(j in 1:J){
  for(k in 1:K){
    n[j,k] <- y[j,k] + u[j,k] + m[j,k]  #total counts
  }
}


### MODEL 1: N-MIXTURE model: ###

cat("
  model {

    #Likelihoods
    for (j in 1:J) {		# cameras
      Nlocal[j] ~ dpois(lambda*area[j])
      for (k in 1:K) {	# frames
        n[j,k] ~ dbin(p, Nlocal[j])
      }#K
    }#J

    #Derived abundance for the entire area
    Ntotal ~ dpois(lambda*Area)

    # Priors
    lambda ~ dgamma(0.01, 0.01)
    p ~ dbeta(1,1)
    
  }#model", fill=TRUE, file="Nmix.txt")

# JAGS data
jags.data <- list(K=K, J=J, n=n, Area=Area_StateSpace, area=area_camera)
Nst <- sum(n)
inits <- function() list(p=runif(1,.6,0.8), Ntotal=Nst) 

# Parameters monitored
params <- c("p", "lambda", "Nlocal")

#MCMC settings
nc <- 3; nt <- 50;  ni <- 200000;  nb <- 50000;  n.adapt <- 20000

# Call JAGS
out <- jags(jags.data, inits, params, "Nmix.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,
            n.adapt = n.adapt, parallel=T)



### MODEL 2: MARKED N-MIXTURE model: ###   
cat("
  model {
           
    # Likelihoods
    # processes at state-space level
    Ntotal ~ dpois(lambda*Area_StateSpace) 
    Mtotal ~ dbin(phi, Ntotal) #Mtotal is known!
    Utotal <- Ntotal-Mtotal

    # Abundance at camera sites
    for (j in 1:J) {		
      M[j] ~ dbin(delta, Mtotal) 
      U[j] ~ dbin(delta, Utotal) 
      N[j] <- M[j]+ U[j]

      # counts by frame
      for (k in 1:K) {	# frames
        m[j,k] ~ dbin(p*theta, M[j])
        u[j,k] ~ dbin(p*theta, U[j])
        y[j,k] ~ dbin(p*(1-theta), N[j])
      }#K occasions
    }#J traps
  
    # Priors
    lambda ~ dgamma(1, 1)       # individuals per unit space
    p ~ dbeta(1,1)              # detection prob given in superpop.
    phi ~ dbeta(1,1)            # prob marked
    theta ~ dbeta(1,1)          # prob marking status is observed if an individual is counted
    delta ~ dbeta(1,1)
  
    ## Derive metrics
    # mean effective area sampled by a camera: if geographically open, effectiveArea will be > area 
    effectiveArea <- delta*Area_StateSpace

    # superpopulation exposed to sampling at a camera
    realized_N_star <- mean(N[1:J])
    expected_N_star <- lambda*Area_StateSpace*delta

  }#model", fill=TRUE, file="MNmix.txt")

# JAGS data
jags.data <- list(K=K, J=J, m=m, u=u, y=y, Area_StateSpace = Area_StateSpace, Mtotal=Mtotal)

#intial values
M_inits <- apply(m, 1, max)+2
U_inits <- apply(u, 1, max)+2

inits <- function() list(p=runif(1,.6,0.8), 
                         lambda=runif(1,.10,0.30), 
                         U=U_inits,
                         M=M_inits,
                         phi=runif(1,.15,0.30), 
                         theta=runif(1,.6,0.95),
                         delta= runif(1,.001,0.01))

# Parameters monitored
params <- c("lambda", "p","theta","phi", "delta", 
            "effectiveArea", "realized_N_star", "expected_N_star" )

# Call JAGS 
out2 <- jags(jags.data, inits, params, "MNmix.txt", 
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, 
             n.adapt = n.adapt, parallel=T) 






