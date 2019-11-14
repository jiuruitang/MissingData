# population: N = 100,000
N = 100000
set.seed(1010)
gamma = 0.85
X <- rbinom(N, 1,p = exp(gamma)/(1+exp(gamma)))
set.seed(112)
pbeta <- rbeta(N,2,5) 
population <- as.data.frame(cbind(X,pbeta))
pop_mean <- mean(population$X) # 0.70324
pop_sd <- sd(population$X)/sqrt(500) # 0.4568321

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
niter <- 10000
burnin <- 5000

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

GAMMA <- rep(NA,niter-burnin)
ALPHA <- matrix(NA,niter-burnin,2)
X_impute <- matrix(NA,niter-burnin,length(X))

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
  cur_mean <- sum(X*W)/sum(W)
  proposed_mean <- sum(proposed_X*W)/sum(W)
  log.r3 = log(pnorm(proposed_mean,mean = pop_mean,sd = pop_sd)) -  log(pnorm(cur_mean,mean = pop_mean,sd = pop_sd))
  if(log(runif(1))< log.r3){
    X <- proposed_X
    X_accepted = X_accepted + 1
  }
  X <- proposed_X
  GAMMA[i-burnin] <- gamma
  ALPHA[i-burnin,] <- alpha
  X_impute[i-burnin,] <- X
}
alpha_accepted/niter
gamma_accepted/niter
X_accepted/niter
plot(ALPHA[1:niter,2],type = "l")
plot(GAMMA,type = "l")
mean(GAMMA)
colMeans(ALPHA)

# get L = 50 data sets
XL <- X_impute[seq(1,5000,100),]

ans_alpha <- matrix(NA,50,2)
for (i in 1:50){
  test_dat <- as.data.frame(cbind(XL[i,],rownames(sub_dat),W,R),stringsAsFactors = FALSE)
  names(test_dat)[1:2] <- c("X1","id")
  test_dat$X <- as.factor(test_dat$X)
  #test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  #test_dat$str <- as.numeric(test_dat$str)
  #test_dat$Y <- as.factor(test_dat$Y)
  test_dat$R <- as.factor(test_dat$R)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- glm(R~X,data = test_dat, family = binomial(link = "probit"))
  ans_alpha[i,] <- coef(m2)
}
colMeans(ans_alpha)


###############################
###############################
#### Try 2 strata
###############################
# population: N = 100,000
N = 50000
N1 <- N*0.7
N2 <- N*0.3
n <- 5000
n1 <- n*0.3
n2 <- n*0.7

set.seed(1010)
gamma1 = 0.85
gamma2 = 0.1
X <- rep(NA,N)
set.seed(112)
X[1:N1] <- rbinom(N1, 1,p = exp(gamma1)/(1+exp(gamma1)))
X[(N1+1):N] <- rbinom(N2, 1,p = exp(gamma2)/(1+exp(gamma2)))
W <- c(rep(N1/n1,N1),rep(N2/n2,N2))

# missing indicator (MNAR)
alpha0 <- -0.8
alpha1 <- -0.5
R_p <- exp(alpha0+alpha1*X)/(1+exp(alpha0+alpha1*X))
set.seed(2130)
Rx <- rbinom(N,1,p = R_p) # 25.47% missing

pop_dat <- as.data.frame(cbind(X,W,Rx))
pop_mean <- mean(pop_dat$X) # 0.57454
pop_sd <- sd(pop_dat$X)/sqrt(n) # 0.00699212

# sample n = 5000 people
# take a sample of 5000 from population
# 1500 from stratum 1, 3500 from statum 2
set.seed(8533)
getSampled <- c(sample(1:N1,n1),sample((N1+1):N2,n2))
sub_dat <- pop_dat[getSampled,]


# helper functions
log_target_alpha <- function(X,R,alpha){
  return(sum(log(exp(alpha[1]+alpha[2]*X)^(1-R))/(1+exp(alpha[1]+alpha[2]*X))))
}

log_target_gamma <- function(X,gamma){
  return(sum(log((exp(gamma)^X)/(1+exp(gamma)))))
}

library(mvtnorm)
niter <- 10000
burnin <- 5000

# Variable names used in sampling
X <- sub_dat$X   # start at the truth
R <- sub_dat$Rx
n_mis = sum(sub_dat$Rx)
W <- sub_dat$W

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

GAMMA <- rep(NA,niter-burnin)
ALPHA <- matrix(NA,niter-burnin,2)
X_impute <- matrix(NA,niter-burnin,length(X))

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
  cur_mean <- sum(X*W)/sum(W)
  proposed_mean <- sum(proposed_X*W)/sum(W)
  log.r3 = log(pnorm(proposed_mean,mean = pop_mean,sd = pop_sd)) -  log(pnorm(cur_mean,mean = pop_mean,sd = pop_sd))
  if(log(runif(1))< log.r3){
    X <- proposed_X
    X_accepted = X_accepted + 1
  }
  X <- proposed_X
  GAMMA[i-burnin] <- gamma
  ALPHA[i-burnin,] <- alpha
  X_impute[i-burnin,] <- X
}
alpha_accepted/niter
gamma_accepted/niter
X_accepted/niter
plot(ALPHA[1:niter,2],type = "l")
plot(GAMMA,type = "l")
mean(GAMMA)
colMeans(ALPHA)

# get L = 50 data sets
XL <- X_impute[seq(1,5000,100),]

ans_alpha <- matrix(NA,50,2)
for (i in 1:50){
  test_dat <- as.data.frame(cbind(XL[i,],rownames(sub_dat),W,R),stringsAsFactors = FALSE)
  names(test_dat)[1:2] <- c("X1","id")
  test_dat$X <- as.factor(test_dat$X)
  #test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  #test_dat$str <- as.numeric(test_dat$str)
  #test_dat$Y <- as.factor(test_dat$Y)
  test_dat$R <- as.factor(test_dat$R)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- glm(R~X,data = test_dat, family = binomial(link = "probit"))
  ans_alpha[i,] <- coef(m2)
}
colMeans(ans_alpha)
