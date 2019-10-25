# population: N = 100,000
N = 100000
set.seed(1010)
gamma = 0.85
X <- rbinom(N, 1,p = exp(gamma)/(1+exp(gamma)))
set.seed(112)
pbeta <- rbeta(N,2,5) 
population <- as.data.frame(cbind(X,pbeta))
pop_mean <- mean(population$X) # 0.70324
pop_sd <- sd(population$X) # 0.4568321

# sample n = 500 people
n <- 500
set.seed(1245)
index_500<- sample(1:N,n,prob = pbeta )
X_500 <- population[index_500,]
X_500$weight <- 1/X_500$pbeta
sum(X_500$weight*X_500$X)/sum(X_500$weight) # 0.67945, quite close to population truth

# missing indicator (MNAR)
alpha0 <- -0.8
alpha1 <- -0.5
R_p <- exp(alpha0+alpha1*X)/(1+exp(alpha0+alpha1*X))
set.seed(2130)
Rx <- rbinom(n,1,p = R_p) # 23.8% missing
# generate observed data
dat_obs <- X_500
dat_obs$Rx <- Rx
dat_obs$X <- ifelse(dat_obs$Rx == 0, dat_obs$X, NA) # Rx = 1 -> X missing
sum(dat_obs$X[which(dat_obs$Rx == 0)]*dat_obs$weight[which(dat_obs$Rx == 0)])/sum(dat_obs$weight[which(dat_obs$Rx == 0)])
#  0.6812081
sum(X_500$X[which(dat_obs$Rx != 0)]*dat_obs$weight[which(dat_obs$Rx != 0)])/sum(dat_obs$weight[which(dat_obs$Rx != 0)])
# 0.6735617

# update other parameters
# when inputing X: 
# do :
#  bootstrap a new data set with prob proportional to weight, with replacement
#  Incorporating the marginal infor into priors, or solve some optimization problem
#  take the majority vote to input a certain missing Xi, if multiple Xi are in the bootstrapped data set
#  calculate the likelihood of getting this imputation,
#  compare to current imputation, decide accept/reject

# helper functions
log_target_alpha <- function(X,R,alpha){
  return(sum(log(exp(alpha[1]+alpha[2]*X)^(1-R))/(1+exp(alpha[1]+alpha[2]*X))))
}

log_target_gamma <- function(X,gamma){
  return(sum(log((exp(gamma)^X)/(1+exp(gamma)))))
}

library(mvtnorm)
niter <- 50000

# Variable names used in sampling
X <- X_500$X   # start at the truth
R <- dat_obs$Rx
n_mis = sum(dat_obs$Rx)
W <- dat_obs$weight

# initial values
alpha = c(-0.8,-0.5)
gamma = 0.85
sigma_alpha0 = 0.25
sigma_alpha1 = 0.25
sigma_gamma = 0.4

# acceptance rate
alpha_accepted = 0
gamma_accepted = 0
X_accepted = 0

GAMMA <- rep(NA,niter)
ALPHA <- matrix(NA,niter,2)
X_impute <- matrix(NA,niter,length(X))

for (i in 1:niter){
  # update gamma
  proposed_gamma = rnorm(1,mean = gamma, sd = sigma_gamma)
  log.r1 = log_target_gamma(X,proposed_gamma) - log_target_gamma(X,gamma)
  if(log(runif(1))< log.r1){
    gamma <- proposed_gamma
    gamma_accepted = gamma_accepted + 1
  }
  
  # update alpha
  # proposed_alpha = rmvnorm(1,mean = current_alpha, sigma = diag(c(sigma_alpha0,sigma_alpha1)))
  proposed_alpha <- rep(NA,2)
  proposed_alpha[1] = rnorm(1,alpha[1],sd = sigma_alpha0)
  proposed_alpha[2] = rnorm(1,alpha[2],sd = sigma_alpha1)
  log.r2 = log_target_alpha(X,R,proposed_alpha) - log_target_alpha(X,R,alpha)
  if(log(runif(1))< log.r2){
    alpha <- proposed_alpha
    alpha_accepted = alpha_accepted + 1
  }
  
  # update X
  # sample an X_mis (R = 1)
  proposed_X <- X
  proposed_X[which(R != 0)] <- rbinom(n_mis, 1,pexp(gamma)*(1+exp(alpha0))/(1+exp(alpha0+alpha1)))
  # M-H based on auxilliary information
  #cur_mean <- sum(X*W)/sum(W)
  #proposed_mean <- sum(proposed_X*W)/sum(W)
  #log.r3 = log(pnorm(proposed_mean,mean = pop_mean,sd = pop_sd)) -  log(pnorm(cur_mean,mean = pop_mean,sd = pop_sd))
  #if(log(runif(1))< log.r3){
  #  X <- proposed_X
  #  X_accepted = X_accepted + 1
  #}
  X <- proposed_X
  GAMMA[i] <- gamma
  ALPHA[i,] <- alpha
  X_impute[i,] <- X
}
alpha_accepted/niter
gamma_accepted/niter
X_accepted/niter
plot(ALPHA[1:niter,2],type = "l")
plot(GAMMA,type = "l")
mean(GAMMA)
colMeans(ALPHA)
