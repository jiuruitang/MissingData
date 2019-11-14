# Implement Michael's stratified sampling

N = 50000
N1 <- N*0.7
N2 <- N*0.3
n <- 5000
n1 <- n*0.3
n2 <- n*0.7

# 1 to 35000 people in stratum 1, 35001 to 50000 people in stratum 2
# theta = (0.5,0.15,0.35) stratum 1, theta = (0.1,0.45,0.45) stratum 2

# sample Y for multinomial
Y <- rep(NA,N)
set.seed(1233)
Y[1:N1] <- sample(c(1,2,3),N1,replace = TRUE, prob = c(0.5,0.15,0.35))
set.seed(2324)
Y[(N1+1):N] <- sample(c(1,2,3),N2,replace = TRUE, prob = c(0.1,0.45,0.45))
W <- c(rep(N1/n1,N1),rep(N2/n2,N2))

# scenario 1
alpha0 <- 0.5
alpha12 <- -0.5
alpha13 <- -1
gamma0 <- -0.25
gamma12 <- 0.1
gamma13 <- 0.3
gamma2 <- -1.1

# sample X1|Y
set.seed(45343)
pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) # pnorm is CDF or normal
set.seed(3204)
X1 <- rbinom(N,1,p=pi_x1)

# sample Rx|X,Y
set.seed(894)
pi_Rx <- pnorm(gamma0 + gamma12*(Y == 2) + gamma13*(Y==3) + gamma2*X)
set.seed(3901)
Rx <- rbinom(N,1,p = pi_Rx)

# get a dataframe for population
pop_dat <- as.data.frame(cbind(Y,X1,Rx,W))

# take a sample of 5000 from population
# 1500 from stratum 1, 3500 from statum 2
set.seed(8533)
getSampled <- c(sample(1:N1,n1),sample((N1+1):N2,n2))
sub_dat <- pop_dat[getSampled,]

# calculate population Total
mean_X1 <- mean(pop_dat$X1)
sd_X1 <- sd(pop_dat$X1)

############################################
############################################
# helper functions for M-H
log_target_alpha <- function(X,Y,alpha){
  #pi_x11 <- log(pnorm(alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3)))
  #pi_x10 <- log(1-pnorm(alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3)))
  #return(sum(pi_x11[which(X==1)]) + sum((1-pi_x10)[which(X==0)]))
  pi_x1 <- pnorm(alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3))
  return(sum(log(pi_x1^X*(1-pi_x1)^(1-X))))
}

log_target_gamma <- function(X,Y,R,gamma){
  #pi_R1 <- log(pnorm(gamma[1] + gamma[2]*(Y==2)+ gamma[3]*(Y==3)+gamma[4]*(X==1)))
  #pi_R0 <- log(1-pnorm(gamma[1] + gamma[2]*(Y==2)+ gamma[3]*(Y==3)+gamma[4]*(X==1)))
  #return(sum(pi_R1[which(R==1)]) + sum((1-pi_R0)[which(R==0)]))
  pi_R1 <- pnorm(gamma[1] + gamma[2]*(Y==2)+ gamma[3]*(Y==3)+gamma[4]*(X==1))
  return(sum(log(pi_R1^R*(1-pi_R1)^(1-R))))
  
}

# Gibbs
X1 <- sub_dat$X1
# X1[which(sub_dat$Rx == 1)] <- NA
W <- sub_dat$W
Rx <- sub_dat$Rx
Y <- sub_dat$Y
n_mis <- sum(sub_dat$Rx == 1)

alpha0 <- 0.5
alpha12 <- -0.5
alpha13 <- -1
gamma0 <- -0.25
gamma12 <- 0.1
gamma13 <- 0.3
gamma2 <- -1.1

alpha <- c(alpha0,alpha12,alpha12)
gamma <- c(gamma0,gamma12,gamma13,gamma2)

sigma_alpha1 = 0.05
sigma_alpha2 = 0.05
sigma_alpha3 = 0.05
sigma_gamma1 =0.05
sigma_gamma2 =0.05
sigma_gamma3 =0.05
sigma_gamma4 =0.05

niter <- 10000
burnin <- 5000

# acceptance ration
X1_accepted = 0
alpha_accepted <- gamma_accepted <- 0

GAMMA <- matrix(NA,niter-burnin,4)
ALPHA <- matrix(NA,niter-burnin,3)
X1_impute <- matrix(NA,niter-burnin,length(X1))

# prior on thetat
theta0 <- c(1,1,1)

for (i in 1:niter){
  # do we need to?
  # update theta
  # theta ~ rdirichlet(1, alpha = c(sum(Y==1),sum(Y==2),sum(Y==3) + theta0))
  
  # for Rx = 1, sample X1_mis*
  pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) 
  pi_Rx_temp <- gamma0 + gamma12*(Y==2) + gamma13*(Y==3)
  q1 <- pi_x1/(1-pi_x1)*pnorm(pi_Rx_temp+gamma2)/pnorm(pi_Rx_temp)
  proposed_X1 <- X1
  proposed_X1[which(Rx==1)] <- rbinom(n_mis, 1, p = q1)
  
  # M-H step
  cur_mean <- sum(X1*W)/sum(W)
  proposed_mean <- sum(proposed_X1*W)/sum(W)
  log.r1 = log(pnorm(proposed_mean,mean = pop_mean,sd = pop_sd)) -  log(pnorm(cur_mean,mean = pop_mean,sd = pop_sd))
  if(log(runif(1))< log.r1){
    X1 <- proposed_X1
    X1_accepted = X1_accepted + 1
  }
  
  # update alpha
  proposed_alpha <- rep(NA,3)
  proposed_alpha[1] = rnorm(1,alpha[1],sd = sigma_alpha1)
  proposed_alpha[2] = rnorm(1,alpha[2],sd = sigma_alpha2)
  proposed_alpha[3] = rnorm(1,alpha[3],sd = sigma_alpha3)
  log.r2 = log_target_alpha(X1,Y,proposed_alpha) - log_target_alpha(X1,Y,alpha)
  if(log(runif(1))< log.r2){
    alpha <- proposed_alpha
    alpha_accepted = alpha_accepted + 1
  }
  
  # update gamma
  proposed_gamma <- rep(NA,4)
  proposed_gamma[1] = rnorm(1,gamma[1],sd = sigma_gamma1)
  proposed_gamma[2] = rnorm(1,gamma[2],sd = sigma_gamma2)
  proposed_gamma[3] = rnorm(1,gamma[3],sd = sigma_gamma3)
  proposed_gamma[4] = rnorm(1,gamma[4],sd = sigma_gamma4)
  log.r3 = log_target_gamma(X1,Y,R,proposed_gamma) - log_target_gamma(X1,Y,R,gamma)
  if(log(runif(1))< log.r3){
    gamma <- proposed_gamma
    gamma_accepted = gamma_accepted + 1
  }
  
  if (i > burnin){
    GAMMA[i-burnin,] <- gamma
    ALPHA[i-burnin,] <- alpha
    X1_impute[i-burnin,] <- X1
  }
}

alpha_accepted/niter
gamma_accepted/niter
X1_accepted/niter  

# get L = 50 data sets
XL <- X1_impute[seq(1,5000,100),]

# design-based estimates
library(survey)
str <- c(rep(1,n1),rep(2,n2))
test_dat <- as.data.frame(cbind(XL[1,],rownames(sub_dat),W,str,Y,Rx),stringsAsFactors = FALSE)
names(test_dat)[1:2] <- c("X1","id")
test_dat$X1 <- as.factor(test_dat$X1)
test_dat$id <- as.numeric(test_dat$id)
test_dat$W <- as.numeric(test_dat$W)
test_dat$str <- as.numeric(test_dat$str)
test_dat$Y <- as.factor(test_dat$Y)
test_dat$Rx <- as.factor(test_dat$Rx)
mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W,strata = ~str)
m1 <- svyglm(X1~Y,design = mydesign,family = binomial(link = "probit"))
summary(m1)

ans <- matrix(NA,50,3)
for (i in 1:50){
  test_dat <- as.data.frame(cbind(XL[i,],rownames(sub_dat),W,str,Y),stringsAsFactors = FALSE)
  names(test_dat)[1:2] <- c("X1","id")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$str <- as.numeric(test_dat$str)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W,strata = ~str)
  m1 <- svyglm(X1~Y,design = mydesign,family = binomial(link = "probit"))
  ans[i,] <- coef(m1)
}
colMeans(ans)

# calculate gamma using stats package
m2 <- glm(Rx~Y+X1,data = test_dat, family = binomial(link = "probit"))
ans_g <- matrix(NA,50,4)
for (i in 1:50){
  test_dat <- as.data.frame(cbind(XL[i,],rownames(sub_dat),W,str,Y,Rx),stringsAsFactors = FALSE)
  names(test_dat)[1:2] <- c("X1","id")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$str <- as.numeric(test_dat$str)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W,strata = ~str)
  m2 <- glm(Rx~Y+X1,data = test_dat, family = binomial(link = "probit"))
  ans_g[i,] <- coef(m2)
}
colMeans(ans_g)

