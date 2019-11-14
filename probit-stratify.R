# Created : 10/29
# Updated: 11/04
#########################################
### Use probit to update alpha and gamma
#########################################
library(tmvtnorm)
library(truncnorm)
library(mvtnorm)

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
#set.seed(1233)
Y[1:N1] <- sample(c(1,2,3),N1,replace = TRUE, prob = c(0.5,0.15,0.35))
#set.seed(2324)
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
#set.seed(45343)
pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) # pnorm is CDF or normal
#set.seed(3204)
X1 <- rbinom(N,1,p=pi_x1)

# sample Rx|X,Y
#set.seed(894)
pi_Rx <- pnorm(gamma0 + gamma12*(Y == 2) + gamma13*(Y==3) + gamma2*X1)
#set.seed(3901)
Rx <- rbinom(N,1,p = pi_Rx)



# get a dataframe for population
pop_dat <- as.data.frame(cbind(Y,X1,Rx,W))
m1 <- glm(Rx~as.factor(Y)+as.factor(X1),data=pop_dat,family=binomial(probit))
summary(m1)

# take a sample of 5000 from population
# 1500 from stratum 1, 3500 from statum 2
#set.seed(8533)
getSampled <- c(sample(1:N1,n1),sample((N1+1):N,n2))
sub_dat <- pop_dat[getSampled,]
m2 <- glm(Rx~as.factor(Y)+as.factor(X1),data=sub_dat,family=binomial(probit))
summary(m2)

# calculate population Total
mean_X1 <- mean(pop_dat$X1)
sd_X1 <- sd(pop_dat$X1)

############################################
############################################
# Gibbs
X1 <- sub_dat$X1
# X1[which(sub_dat$Rx == 1)] <- NA
W <- sub_dat$W
Rx <- sub_dat$Rx
Y <- sub_dat$Y
n_mis <- sum(sub_dat$Rx == 1)
m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)

# auxilliary info
pop_mean <- mean(pop_dat$X1)
pop_sd <- sd(pop_dat$X1)/sqrt(n)


alpha <- rnorm(3)
gamma <- rnorm(4)

niter <- 10000
burnin <- 5000

# acceptance ration
X1_accepted = 0

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
  proposed_X1[which(Rx==1)] <- rbinom(n_mis, 1, p = q1/(1+q1))
  
  # M-H step
  cur_mean <- sum(X1*W)/sum(W)
  proposed_mean <- sum(proposed_X1*W)/sum(W)
  log.r1 = log(dnorm(proposed_mean,mean = pop_mean,sd = pop_sd)) -  log(dnorm(cur_mean,mean = pop_mean,sd = pop_sd))
  if(log(runif(1))< log.r1){
    X1 <- proposed_X1
    X1_accepted = X1_accepted + 1
  }
  
  # update alpha, with data augumentation
  # update Z|X,alpha
  Z <- X1
  Y_n <- Y[which(X1 == 0)]
  Z_mean_n <- alpha[1] + alpha[2]*(Y_n==2) + alpha[3]*(Y_n==3)
  Z[which(X1 == 0)] <- rtruncnorm(1, a=-Inf, b=0, mean = Z_mean_n, sd = 1)
  Y_p <- Y[which(X1 != 0)]
  Z_mean_p <- alpha[1] + alpha[2]*(Y_p==2) + alpha[3]*(Y_p==3)
  Z[which(X1 != 0)] <- rtruncnorm(1, a=0, b=Inf, mean = Z_mean_p, sd = 1)
  # update alpha|Z,Y
  Y_mat <- model.matrix(X1~as.factor(Y))
  sigma_hat <- solve(diag(length(alpha)) + t(Y_mat)%*%Y_mat)
  alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z)
  alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
  
  # update gamma
  # update G|Rx, gamma
  G <- Rx
  Y_n <- Y[which(Rx == 0)]
  X_n <- X1[which(Rx == 0)]
  G_mean_n <- gamma[1] + gamma[2]*(Y_n==2) + gamma[3]*(Y_n==3) + gamma[4]*X_n
  G[which(Rx == 0)] <- rtruncnorm(1, a=-Inf, b=0, mean = G_mean_n, sd = 1)
  Y_p <- Y[which(Rx != 0)]
  X_p <- X1[which(Rx != 0)]
  G_mean_p <- gamma[1] + gamma[2]*(Y_p==2) + gamma[3]*(Y_p=3) + gamma[4]*X_p
  G[which(Rx != 0)] <- rtruncnorm(1, a=0, b=Inf, mean = G_mean_p, sd = 1)
  # update alpha|Z,Y
  YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
  sigma_hat <- solve(diag(length(gamma)) + t(YX_mat)%*%YX_mat)
  gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G)
  gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
  
  
  if (i > burnin){
    GAMMA[i-burnin,] <- gamma
    ALPHA[i-burnin,] <- alpha
    X1_impute[i-burnin,] <- X1
  }
}


X1_accepted/niter  

# get L = 50 data sets
XL <- X1_impute[seq(1,5000,100),]

# design-based estimates
library(survey)

ans <- matrix(NA,50,3)
for (i in 1:50){
  str <- c(rep(1,n1),rep(2,n2))
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
# 0.5,-0.5,-1

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
# -0.25, 0.1, 0.3,-1.1

################################
############### 11/04 Update
library(tmvtnorm)
library(truncnorm)
library(mvtnorm)

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
#set.seed(1233)
Y[1:N1] <- sample(c(1,2,3),N1,replace = TRUE, prob = c(0.5,0.15,0.35))
#set.seed(2324)
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
#set.seed(45343)
pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) # pnorm is CDF or normal
#set.seed(3204)
X1 <- rbinom(N,1,p=pi_x1)

# sample Rx|X,Y
#set.seed(894)
pi_Rx <- pnorm(gamma0 + gamma12*(Y == 2) + gamma13*(Y==3) + gamma2*X1)
#set.seed(3901)
Rx <- rbinom(N,1,p = pi_Rx)



# get a dataframe for population
pop_dat <- as.data.frame(cbind(Y,X1,Rx,W))
m1 <- glm(Rx~as.factor(Y)+as.factor(X1),data=pop_dat,family=binomial(probit))
summary(m1)

# take a sample of 5000 from population
# 1500 from stratum 1, 3500 from statum 2
#set.seed(8533)
getSampled <- c(sample(1:N1,n1),sample((N1+1):N,n2))
sub_dat <- pop_dat[getSampled,]
m2 <- glm(Rx~as.factor(Y)+as.factor(X1),data=sub_dat,family=binomial(probit))
summary(m2)

# calculate population Total
mean_X1 <- mean(pop_dat$X1)
sd_X1 <- sd(pop_dat$X1)

############################################
############################################
# Gibbs
X1 <- sub_dat$X1
# X1[which(sub_dat$Rx == 1)] <- NA
W <- sub_dat$W
Rx <- sub_dat$Rx
Y <- sub_dat$Y
n_mis <- sum(sub_dat$Rx == 1)
m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)

# auxilliary info
pop_mean <- mean(pop_dat$X1)
pop_sd <- sd(pop_dat$X1)/sqrt(n)

sub_dat$Y <- as.factor(sub_dat$Y)
sub_dat$X1 <- as.factor(sub_dat$X1)

Model <- list(X1=~Y,Rx=~Y+X1)
alpha <- rnorm(ncol(model.matrix(Model$X1,data=sub_dat)))
gamma <- rnorm(ncol(model.matrix(Model$R,data=sub_dat)))
b0_alpha <- alpha; b0_alpha[] <- 0
sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
b0_gamma <- gamma; b0_gamma[] <- 0
sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
Z <- data.frame(X1=rnorm(nrow(sub_dat)))
Z$R <- rnorm(nrow(sub_dat))
  
niter <- 10000
burnin <- 5000

# acceptance ration
X1_accepted = 0

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
  proposed_X1[which(Rx==1)] <- rbinom(n_mis, 1, p = q1/(1+q1))
  
  # M-H step
  cur_mean <- sum(X1*W)/sum(W)
  proposed_mean <- sum(proposed_X1*W)/sum(W)
  log.r1 = log(dnorm(proposed_mean,mean = pop_mean,sd = pop_sd)) -  log(dnorm(cur_mean,mean = pop_mean,sd = pop_sd))
  if(log(runif(1))< log.r1){
    X1 <- proposed_X1
    X1_accepted = X1_accepted + 1
  }
  
  ## sample alpha, the parameters for X1|Y1
  X1_cond <- model.matrix(Model$X1,data=sub_dat)
  mu_alpha <- solve((t(X1_cond)%*%X1_cond)+solve(sigma0_alpha))%*%
    ((t(X1_cond)%*%Z$X1)+(solve(sigma0_alpha)%*%b0_alpha))
  sigma_alpha <- solve((t(X1_cond)%*%X1_cond)+solve(sigma0_alpha))
  alpha <- rmvnorm(1,mean=mu_alpha,sigma=sigma_alpha)
  
  ## sample Z_X1, the augmented variables for X1
  mu_Z_X1 <- X1_cond%*%t(alpha)
  U_X1 <- mu_Z_X1; U_X1[] <- 0
  U_X1[sub_dat$X1==0] <- runif(sum(sub_dat$X1==0),pnorm(-Inf-mu_Z_X1[sub_dat$X1==0]),pnorm(0-mu_Z_X1[sub_dat$X1==0]))
  U_X1[sub_dat$X1==1] <- runif(sum(sub_dat$X1==1),pnorm(0-mu_Z_X1[sub_dat$X1==1]),pnorm(Inf-mu_Z_X1[sub_dat$X1==1]))
  Z$X1 <- mu_Z_X1 + qnorm(U_X1)
  
  ## sample gamma, the parameters for R|Y1,X1
  R_cond <- model.matrix(Model$R,data=sub_dat)
  mu_gamma <- solve((t(R_cond)%*%R_cond)+solve(sigma0_gamma))%*%
    ((t(R_cond)%*%Z$R)+(solve(sigma0_gamma)%*%b0_gamma))
  sigma_gamma <- solve((t(R_cond)%*%R_cond)+solve(sigma0_gamma))
  gamma <- rmvnorm(1,mean=mu_gamma,sigma=sigma_gamma)
  
  ## sample Z_R, the augmented variables for R
  mu_Z_R <- R_cond%*%t(gamma)
  U_Z <- mu_Z_R; U_Z[] <- 0
  U_Z[sub_dat$Rx==0] <- runif(sum(sub_dat$Rx==0),pnorm(-Inf-mu_Z_R[sub_dat$Rx==0]),pnorm(0-mu_Z_R[sub_dat$Rx==0]))
  U_Z[sub_dat$Rx==1] <- runif(sum(sub_dat$Rx==1),pnorm(0-mu_Z_R[sub_dat$Rx==1]),pnorm(Inf-mu_Z_R[sub_dat$Rx==1]))
  Z$R <- mu_Z_R + qnorm(U_Z)
  
  
  if (i > burnin){
    GAMMA[i-burnin,] <- gamma
    ALPHA[i-burnin,] <- alpha
    X1_impute[i-burnin,] <- X1
  }
}


X1_accepted/niter  

# get L = 50 data sets
XL <- X1_impute[seq(1,5000,100),]

# design-based estimates
library(survey)

ans <- matrix(NA,50,3)
for (i in 1:50){
  str <- c(rep(1,n1),rep(2,n2))
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
# 0.5,-0.5,-1

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
# -0.25, 0.1, 0.3,-1.1