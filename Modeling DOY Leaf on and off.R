### Leaf on (Lon) and leaf off (Loff) models

# load required packages
library(jagsUI) # must install jags program to your computer first
library(MCMCpack)
library(ggmcmc)
library(coda)
library(beepr)
library(data.table)

##################################
##################################
# Lon and latitude model... 

# Set working directory to where the data are saved/stored
setwd("insert wd here")

# Load the data
Lon <- fread("Lon.csv", na.strings="NA",header=T)
colnames(Lon)

# Standardize covariates in new matrix called dat (for latitude and env. var. models)
dat <- Lon
changeCols <- colnames(dat)[c(3,10:12)]
dat[,(changeCols):= lapply(.SD, scale), .SDcols = changeCols]
summary(dat)

# Covariate matrix - select desired covariate
x_mat <- as.matrix(dat[,c("Latitude")])
head(x_mat)

################################################################
## MODEL for Lon or Loff by Latitude and species nested within nativity with year as a random variable

# Define the model in the jags language and write a text file
sink("1coefficient.txt")
cat("
    model {
    
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dnorm(mu[i], tau[group[i]])               
    mu[i] <- alpha[group[i]] + beta[group[i]] * x[i,1] + delta[year[i]] 
    e.y[i] <- y[i] - mu[i]
    } 
    
    # Hierarchical variance model 
    for(j in 1:J) { # number of species
    log.sigma[j] ~ dnorm(mu.sigma[j], inv.omega.sigma.squared)
    log(sigma[j]) <- log.sigma[j]
    tau[j] <- 1/pow(sigma[j], 2)
    mu.sigma[j] <- Mu.Sigma  + gamma.sigma[1] * z[j]
    }
    
    # Priors for variances
    Mu.Sigma ~ dnorm(0, 0.001)
    Med.Sigma <- exp(Mu.Sigma)
    omega.sigma ~ dunif(0, 100)
    inv.omega.sigma.squared <- 1/pow(omega.sigma, 2)
    
    # Normal priors on slopes for modeling the residual SD
    for(i in 1:1){
    gamma.sigma[i] ~ dnorm(0, 0.0001)
    }
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivriate normal
    
    BB.hat[j,1] <- mu.alpha + gamma1a * z[j] 
    e.a[j] <- alpha[j] - BB.hat[j,1]
    
    BB.hat[j,2] <- mu.beta[1]  + gamma1b[1] * z[j] 
    e.b[j] <- beta[j] - BB.hat[j,2]
    
    }
    
    # Priors 
    for(i in 1:ncovs){
    mu.beta[i] ~ dnorm(0, 0.001)
    gamma1b[i] ~ dnorm(0, 0.001)
    }
    
    # Slope parameters for how (native and exotic) modify the covariate-DOY relationship
    mu.alpha ~ dnorm(0, 0.0001) 
    gamma1a ~ dnorm(0, 0.0001)    
    
    for(m in 1:M){
    delta[m] ~ dnorm(0, tau.m)
    }
    
    sigma.m ~ dunif(0,100)
    tau.m <- pow(sigma.m,-2)
    sigma.m2 <- pow(sigma.m,2)    
    
    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()

# Number of species
J <- length(unique(dat$Common_Name))

# Species indicator
dat$G <- as.numeric(as.factor(as.numeric(as.factor(dat$Common_Name))))

# Make nativity a binary variable 
dat$nat_binary <- ifelse(dat$Nativity == 'exotic', 1, 0)
nativity <- as.numeric(by(dat$nat_binary, dat$Common_Name, mean))

# Converts year into a factor and then back to numeric, making it a number factor
dat$year <- as.numeric(as.factor(as.numeric(as.factor(dat$First_Yes_Year))))   
# Assigns year to each observation
m <- as.numeric(by(dat$year, dat$First_Yes_Year, mean))                        
# Object holding the number of years
M = length(unique(dat$year))                                          

# Number of varying parameters
K <- 2

# Create identity matrix for Wishart dist'n
W <- diag(K)

# Load data
data <- list(y = dat$First_Yes_DOY, group = dat$G, n = dim(dat)[1], J = J,
             x=x_mat, K=K, W=W, z = nativity, M=M, year = dat$year,ncovs = dim(x_mat)[2])

# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1),  
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma1a = rnorm(1), sigma.m = runif(1),
        Mu.Sigma = rnorm(1), omega.sigma = runif(1,0,100))
}

# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B", "gamma1a", "gamma1b", "sigma.m", "delta","sigma",
                "Mu.Sigma", "omega.sigma","gamma.sigma","e.a","e.b",
                "e.y","alpha","beta") 

# MCMC settings
ni <- 120000
nt <- 3
nb <- 90000
nc <- 3

# Call JAGS from R (compute)
Lon_Lat <- jags(data, inits, parameters, "1coefficient.txt", n.chains = nc, 
                n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, DIC = T)
beep(sound=2) # To indicate model is done running

# Look whether model converged (and DIC value)
Lon_Lat
# Save the summary file of the jags output, if desired
# write.csv(Lon_Lat$summary,'Lon_Lat.csv')
# # Save the whole thing
# save(Lon_Lat, file = "Lon_Lat")
# # Load it later
# load("Lon_Lat")

# Calculate R2 
y = dat$First_Yes_DOY
# data level summaries
rsquared.y <- 1 - mean (apply (Lon_Lat$sims.list$e.y, 1, var))/ var (y)

## Look at some model output
# Significance of nativity covariates on varying slopes (relationship with Latitude)
# i.e. Difference in slopes between native and exotic species groups
mean(Lon_Lat$sims.list$gamma1b)
quantile(Lon_Lat$sims.list$gamma1b, c(0.025,0.975))

# Slope for exotics... 
mean(Lon_Lat$sims.list$mu.beta + Lon_Lat$sims.list$gamma1b)
quantile(Lon_Lat$sims.list$mu.beta + Lon_Lat$sims.list$gamma1b, c(0.025,0.975))

# Slope for natives... 
mean(Lon_Lat$sims.list$mu.beta)
quantile(Lon_Lat$sims.list$mu.beta, c(0.025,0.975))

# Test for differences in slope
exotic.slope <- Lon_Lat$sims.list$mu.beta + Lon_Lat$sims.list$gamma1b
native.slope <- Lon_Lat$sims.list$mu.beta
diff <- exotic.slope - native.slope
quantile(diff, c(0.025,0.975)) # 95% CI on difference, does not overlap with zero, so "sig. diff"

##################################
##################################
# Leaf on and environmental variables model... 

# Covariate matrix for the 3 environmental variables for Leaf on
x_mat <- as.matrix(dat[,c("w.precip",	"AGDD",	"chill.days")])
head(x_mat)

# Check for multicollineearity of predictor variables
library(HH)
x_vif <- as.data.frame(x_mat)
vif(x_vif)

################################################################
## MODEL for Lon or Loff with 3 environmental variable and species nested within nativity with year as a random variable

# Define the model in the jags language and write a text file
sink("3coefficients.txt")
cat("
    
    model {
    
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dnorm(mu[i], tau[group[i]])               
    mu[i] <- alpha[group[i]] + beta[group[i]] * x[i,1] + beta2[group[i]] * x[i,2] + beta3[group[i]] * x[i,3] + delta[year[i]] 
    e.y[i] <- y[i] - mu[i]
    } 
    
    # Hierarchical variance model 
    for(j in 1:J) { # number of species
    log.sigma[j] ~ dnorm(mu.sigma[j], inv.omega.sigma.squared)
    log(sigma[j]) <- log.sigma[j]
    tau[j] <- 1/pow(sigma[j], 2)
    mu.sigma[j] <- Mu.Sigma  + gamma.sigma[1] * z[j]
    }
    
    # Priors for variances
    Mu.Sigma ~ dnorm(0, 0.001)
    Med.Sigma <- exp(Mu.Sigma)
    omega.sigma ~ dunif(0, 100)
    inv.omega.sigma.squared <- 1/pow(omega.sigma, 2)
    
    # Normal priors on slopes for modeling the residual SD
    for(i in 1:1){
    gamma.sigma[i] ~ dnorm(0, 0.0001)
    }

    # Level-2 of the model
    for(j in 1:J){ #j is 14 species, group-level predictor
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
    beta2[j] <- BB[j,3]
    beta3[j] <- BB[j,4]
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivriate normal
    
    BB.hat[j,1] <- mu.alpha + gamma1a * z[j] 
    e.a[j] <- alpha[j] - BB.hat[j,1]
    
    BB.hat[j,2] <- mu.beta[1]  + gamma1b[1] * z[j] 
    e.b[j] <- beta[j] - BB.hat[j,2]
    
    BB.hat[j,3] <- mu.beta[2]  + gamma1b[2] * z[j] 
    e.b2[j] <- beta2[j] - BB.hat[j,3]
    
    BB.hat[j,4] <- mu.beta[3]  + gamma1b[3] * z[j] 
    e.b3[j] <- beta3[j] - BB.hat[j,4]
    }
    
    # Priors 
    for(i in 1:ncovs){
    mu.beta[i] ~ dnorm(0, 0.001)
    gamma1b[i] ~ dnorm(0, 0.001)
    }
    
    # Slope parameters for how (native and exotic) modify the covariate-DOY relationship
    mu.alpha ~ dnorm(0, 0.0001) 
    gamma1a ~ dnorm(0, 0.0001)    
    
    for(m in 1:M){
    delta[m] ~ dnorm(0, tau.m)
    }
    
    sigma.m ~ dunif(0,10)
    tau.m <- pow(sigma.m,-2)
    sigma.m2 <- pow(sigma.m,2)    
    
    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()

# Number of species
J <- length(unique(dat$Common_Name))

# Species indicator
dat$G <- as.numeric(as.factor(as.numeric(as.factor(dat$Common_Name))))

# Make nativity a binary variable 
dat$nat_binary <- ifelse(dat$Nativity == 'exotic', 1, 0)
nativity <- as.numeric(by(dat$nat_binary, dat$Common_Name, mean))

# Converts year into a factor and then back to numeric, making it a number factor
dat$year <- as.numeric(as.factor(as.numeric(as.factor(dat$First_Yes_Year))))   
# Assigns year to each observation
m <- as.numeric(by(dat$year, dat$First_Yes_Year, mean))                        
# Object holding the number of years
M = length(unique(dat$year))                                          

# Number of varying parameters
K <- 4
# Create identity matrix for Wishart dist'n
W <- diag(K)

# Load data
data <- list(y = dat$First_Yes_DOY, group = dat$G, n = dim(dat)[1], J = J,
             x = x_mat, K = K, W = W, z = nativity, M = M, year = dat$year, ncovs = dim(x_mat)[2])

# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1),  
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma1a = rnorm(1), sigma.m = runif(1),
        Mu.Sigma = rnorm(1), omega.sigma = runif(1,0,100))
}

# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B", "gamma1a", "gamma1b", "sigma.m", "delta","sigma",
                "Mu.Sigma", "omega.sigma","gamma.sigma","e.a","e.b","e.b2","e.b3",
                "e.y","alpha","beta","beta2","beta3") 

# Call JAGS from R (compute)
Lon_env.var <- jags(data, inits, parameters, "3coefficients.txt", n.chains = nc, 
                                    n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, DIC = T)
beep(sound=8)

# Look whether model converged (and DIC value)
Lon_env.var

# Save the summary file of the jags output, if desired
# write.csv(Lon_env.var$summary,'Lon_env.var.csv')
# # Save the whole thing
# save(Lon_env.var, file = "Lon_env.var")
# # Load it later
# load("Lon_env.var")

##################################
##################################
# Leaf off and latitude model... 

# Set working directory to where the data are saved/stored
setwd("insert wd here")

# Load the data
Loff <- fread("Loff.csv", na.strings="NA",header=T)
colnames(Loff)

# Standardize covariates in new matrix called dat
dat <- Loff
changeCols <- colnames(dat)[c(3,10:14)]
dat[,(changeCols):= lapply(.SD, scale), .SDcols = changeCols]
summary(dat)

# Select latitude covariate
x_mat <- as.matrix(dat[,c("Latitude")])
head(x_mat)

# Number of species
J <- length(unique(dat$Common_Name))

# Species indicator
dat$G <- as.numeric(as.factor(as.numeric(as.factor(dat$Common_Name))))

# Make nativity a binary variable 
dat$nat_binary <- ifelse(dat$Nativity == 'exotic', 1, 0)
nativity <- as.numeric(by(dat$nat_binary, dat$Common_Name, mean))

dat$year <- as.numeric(as.factor(as.numeric(as.factor(dat$Last_Yes_Year))))   #Converts year into a factor and then back to numeric, making it a number factor
m <- as.numeric(by(dat$year, dat$Last_Yes_Year, mean))                        #Assigns year to each observation
M = length(unique(dat$year))                                          #Object holding the number of sites

# Number of varying parameters
K <- 2

# Create identity matrix for Wishart dist'n
W <- diag(K)

# Load data
data <- list(y = dat$Last_Yes_DOY, group = dat$G, n = dim(dat)[1], J = J,
             x = x_mat, K = K, W = W, z = nativity, M = M, year = dat$year, ncovs = dim(x_mat)[2])

# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1),  
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma1a = rnorm(1), sigma.m = runif(1),
        Mu.Sigma = rnorm(1), omega.sigma = runif(1,0,100))
}

# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B", "gamma1a", "gamma1b", "sigma.m", "delta","sigma",
                "Mu.Sigma", "omega.sigma","gamma.sigma","e.a","e.b",
                "e.y","alpha","beta") 

# Call JAGS from R (compute)
Loff_Lat <- jags(data, inits, parameters, "1coefficient.txt", n.chains = nc, 
                 n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, DIC = T)
beep(sound=2)

# Check model convergence, etc... 
Lon_Lat
# Save the summary file of the jags output, if desired
# write.csv(Loff_Lat$summary,'Loff_Lat.csv')
# # Save the whole thing
# save(Loff_Lat, file = "Loff_Lat")
# # Load it later
# load("Loff_Lat")

##################################
##################################
# Leaf off and environmental variables model with temperature separated by season (3 variables)

# Covariate matrix of selected environmental variables for the 3-variable model that separates seasons
x_mat <- as.matrix(dat[,c("summer.P",	"summer.Tmax",	"fall.Tmin")])
head(x_mat)

# Check for multicollineearity of predictor variables
x_vif <- as.data.frame(x_mat)
vif(x_vif)

# Number of species
J <- length(unique(dat$Common_Name))

# Species indicator
dat$G <- as.numeric(as.factor(as.numeric(as.factor(dat$Common_Name))))

# Make nativity a binary variable 
dat$nat_binary <- ifelse(dat$Nativity == 'exotic', 1, 0)
nativity <- as.numeric(by(dat$nat_binary, dat$Common_Name, mean))

# Converts year into a factor and then back to numeric, making it a number factor
dat$year <- as.numeric(as.factor(as.numeric(as.factor(dat$Last_Yes_Year))))   
# Assigns year to each observation
m <- as.numeric(by(dat$year, dat$Last_Yes_Year, mean))                        
# Object holding the number of years
M = length(unique(dat$year))                                          

# Number of varying parameters
K <- 4
# Create identity matrix for Wishart dist'n
W <- diag(K)

# Load data
data <- list(y = dat$Last_Yes_DOY, group = dat$G, n = dim(dat)[1], J = J,
             x = x_mat, K = K, W = W, z = nativity, M = M, year = dat$year, ncovs = dim(x_mat)[2])

# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1),  
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma1a = rnorm(1), sigma.m = runif(1),
        Mu.Sigma = rnorm(1), omega.sigma = runif(1,0,100))
}

# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B", "gamma1a", "gamma1b", "sigma.m", "delta","sigma",
                "Mu.Sigma", "omega.sigma","gamma.sigma","e.a","e.b","e.b2","e.b3",
                "e.y","alpha","beta","beta2","beta3") 

# Call JAGS from R (compute)
Loff_3env.var <- jags(data, inits, parameters, "3coefficients.txt", n.chains = nc, 
                               n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, DIC = T)
beep(sound=2)

# Check for convergence, and look at values
Loff_3env.var
# Save the summary file of the jags output, if desired
# write.csv(Loff_3env.var$summary,'Loff_3env.var.csv')
# # Save the whole thing
# save(Loff_3env.var, file = "Loff_3env.var")
# # Load it later
# load("Loff_3env.var")

##################################
##################################
# Leaf off and environmental variables model across growing season (2 variables)

# Covariate matrix of selected environmental variables for the 3-variable model that separates seasons
x_mat <- as.matrix(dat[,c("GS.P",	"GS.T")])
head(x_mat)

# Check for multicollineearity of predictor variables
x_vif <- as.data.frame(x_mat)
vif(x_vif)

################################################################
## MODEL for 2 environmental variable and species nested within nativity with year as a random variable

# Define the model in the jags language and write a text file
sink("2coefficients.txt")
cat("
    model {
    
    # Level-1 of the model
    for (i in 1:n){ 
       y[i] ~ dnorm(mu[i], tau[group[i]])               
       mu[i] <- alpha[group[i]] + beta[group[i]] * x[i,1] + beta2[group[i]] * x[i,2] + delta[year[i]] 
      e.y[i] <- y[i] - mu[i]
    } 
    
    # Hierarchical variance model 
    for(j in 1:J) { # Number of species
      log.sigma[j] ~ dnorm(mu.sigma[j], inv.omega.sigma.squared)
      log(sigma[j]) <- log.sigma[j]
      tau[j] <- 1/pow(sigma[j], 2)
      mu.sigma[j] <- Mu.Sigma  + gamma.sigma[1] * z[j]
    }
    
    # Priors for variances
    Mu.Sigma ~ dnorm(0, 0.001)
    Med.Sigma <- exp(Mu.Sigma)
    omega.sigma ~ dunif(0, 100)
    inv.omega.sigma.squared <- 1/pow(omega.sigma, 2)

    # Normal priors on slopes for modeling the residual SD
    for(i in 1:1){
      gamma.sigma[i] ~ dnorm(0, 0.0001)
    }


    # Level-2 of the model 
    for(j in 1:J){ 
      alpha[j] <- BB[j,1]
      beta[j] <- BB[j,2]
      beta2[j] <- BB[j,3]
    
      BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivriate normal
    
      BB.hat[j,1] <- mu.alpha + gamma1a * z[j] 
      e.a[j] <- alpha[j] - BB.hat[j,1]

      BB.hat[j,2] <- mu.beta[1]  + gamma1b[1] * z[j] 
      e.b[j] <- beta[j] - BB.hat[j,2]

      BB.hat[j,3] <- mu.beta[2]  + gamma1b[2] * z[j] 
      e.b2[j] <- beta2[j] - BB.hat[j,3]

     }
    
  # Priors 
  for(i in 1:ncovs){
    mu.beta[i] ~ dnorm(0, 0.001)
    gamma1b[i] ~ dnorm(0, 0.001)
  }

  # Slope parameters for how (native and exotic) modify the covariate-DOY relationship
  mu.alpha ~ dnorm(0, 0.0001) 
  gamma1a ~ dnorm(0, 0.0001)    

    for(m in 1:M){
      delta[m]~dnorm(0, tau.m) # This is zero becuase we're not adding another intercept
    }

  sigma.m ~ dunif(0,10)
  tau.m <- pow(sigma.m,-2)
  sigma.m2 <- pow(sigma.m,2)    

    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()

# Number of species
J <- length(unique(dat$Common_Name))

# Species indicator
dat$G <- as.numeric(as.factor(as.numeric(as.factor(dat$Common_Name))))

# Make developmental mode a binary variable 
dat$nat_binary <- ifelse(dat$Nativity == 'exotic', 1, 0)
nativity <- as.numeric(by(dat$nat_binary, dat$Common_Name, mean))

# Converts year into a factor and then back to numeric, making it a number factor
dat$year <- as.numeric(as.factor(as.numeric(as.factor(dat$Last_Yes_Year))))   
# Assigns year to each observation
m <- as.numeric(by(dat$year, dat$Last_Yes_Year, mean))                        
# Object holding the number of years
M = length(unique(dat$year))                                          

# Number of varying parameters
K <- 3
# Create identity matrix for Wishart dist'n
W <- diag(K)

# Load data
data <- list(y = dat$Last_Yes_DOY, group = dat$G, n = dim(dat)[1], J = J,
             x=x_mat, K=K, W=W, z = nativity, M=M, year = dat$year,ncovs = dim(x_mat)[2])

# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1),  
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),Tau.B=rwish(K+1,diag(K)),
        gamma1a = rnorm(1), sigma.m = runif(1),
        Mu.Sigma = rnorm(1), omega.sigma = runif(1,0,100))
}

# Parameters monitored
parameters <- c("mu.alpha","mu.beta","BB", "Sigma.B", "gamma1a", "gamma1b", "sigma.m", "delta","sigma",
                "Mu.Sigma", "omega.sigma","gamma.sigma","e.a","e.b","e.b2",
                "e.y","alpha","beta","beta2") 

# Call JAGS from R (compute)
Loff_2env.var <- jags(data, inits, parameters, "2coefficients.txt", n.chains = nc, 
                      n.thin = nt, n.iter = ni, n.burnin = nb, parallel=T, DIC = T)
beep(sound=2)

# Check for convergence, and look at values
Loff_2env.var
# Save the summary file of the jags output, if desired
# write.csv(Loff_2env.var$summary,'Loff_2env.var.csv')
# # Save the whole thing
# save(Loff_2env.var, file = "Loff_2env.var")
# # Load it later
# load("Loff_2env.var")
