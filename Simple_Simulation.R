# Created 09/26
# Modified 09/26
# This code aims to generate a simulation data set 
# for missing data with auxiliary infor
# and survey weights

# 500 people
n <- 500

# sample Y for bernoulli p = 0.7
Y <- rbinom(500,1,prob = 0.7)

# sample X|Y fro bernoulli pi
# where pi is normal(alpha0+alpha1 + 1(Y1=1))
alpha0 = 0.5
alpha1 = -0.5
pi <- pnorm(alpha0 + alpha1*Y)

X <- rbinom(500,1,prob = pi)

# sample weight from lognormal
Z <- rlnorm(500, meanlog = 0, sdlog = 1)
pi_survey <- Z*n/sum(Z) # prob of selection i-th unit
pi_survey <- ifelse(pi_survey > 1, 1,pi_survey) # make prob equals 1 if greater than 1
W <- 1/pi_survey # weight is one over prob

# auxiliary information
sum(W*X) # 817

# generate missing pattern
gamma0 <- -0.4
gamma1 <- 0.1
gamma2 <- -1.5

piR <- pnorm(gamma0 + gamma1*Y + gamma2*X)
R <- rbinom(n,1,prob = piR) #18 ~ 20% missing

sum(W[which(R==0)]*X[which(R==0)]) # 796
