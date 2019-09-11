# Created 09/07
# This R script aim to replicate the AN model in Deng's master thesis
# In chapter2, with MCMC algorithm in Appendix A

###############################################
# Defind the parameters in the generating model
beta0 <- 0.3
beta1 <- -0.4
gamma0 <- 0.3
gamma1 <- -0.3
gamma2 <- 0.7
alpha0 <- -0.4
alpha1 <- 1.0
alpha2 <- -0.7
alpha3 <- 1.3

# my choose for total sample size N, and p
## total number  = Ncp + Nip + Nr = 1100, observe X
## Np = Ncp + Nip = 800, observe X,Y1
## Ncp + Nr, observe Y2
p <- 0.7
Np <- 800
Nr <- 300


# for i = 1,...,Np, observe Xi and Y1i
X <- rbinom(n = Np + Nr, size = 1, prob = p)

Y1_p <- exp(beta0 + beta1*X)/(1+exp(beta0+beta1*X))
Y1 <- rbinom(n = Np + Nr, size = 1, prob = Y1_p)

Y2_p <- exp(gamma0+gamma1*X+gamma2*Y1)/(1+exp(gamma0+gamma1*X+gamma2*Y1))
Y2 <- rbinom(n = Np + Nr, size  = 1, prob = Y2_p)

W1_p <- exp(alpha0+alpha1*X+alpha2*Y1+alpha3*Y2)/(1+exp(alpha0+alpha1*X+alpha2*Y1+alpha3*Y2))
W1 <- rbinom(n = Np + Nr, size  = 1, prob = W1_p)

# for Y1, only observe Ncp + Nip of them
Y1_obs <- c(Y1[1:(Np)],rep(NA,Nr))

# for W1, only observe Ncp + Nip of them
W1_obs <- c(W1[1:Np],rep(NA,Nr))

# for Y2, missing when W1 = 0
Y2_obs <- Y2
Y2_obs[which(W1_obs == 0)] <- NA


##############################################
##### MCMC algorithm

# supportive functions
iter <- 10000
burnIn <- 5000

# starting value
g0 <- g1 <- g2 <- b0 <- b1 <- a0 <- a1 <- a2 <- a3 <- 0.1

for (i in 1:iter){
  # Step 1
  # for units with W1i = 0, impute Y2i, based on observed Xi and Y1i
  q1 <- exp(g0 + g1*X + g2*Y1_obs)*(1+exp(a0+a1*X+a2*Y1_obs))/(1+exp(a0+a1*X+a2*Y1_obs+a3))
  q_step1 <- (q1/(1+q1))[which(W1_obs == 0)]
  Y2_obs[which(W1_obs == 0)] <- rbinom(n = length(which(W1_obs == 0)), size  = 1, prob = q_step1)
  
  # Step 2
  # for units in the refreshment sample, impute Y1i, based on X, Y2 and W1
  q2 <- exp(b0+b1*X+g2*Y2_obs+a2*W1_obs)*(1+exp(g0+g1*X))/(1+exp(g0+g1*X+g2))*(1+exp(a0+a1*X+a3*Y2_obs))/(1+exp(a0+a1*X+a2+a3*Y2_obs))
  q_step2 <- (q2/(1+q2))[(Np+1):(Np+Nr)]
  n_step2 <- length(Y1_obs[(Np+1):(Np+Nr)])
  Y1_obs[(Np+1):(Np+Nr)] <- rbinom(n = n_step2, size = 1, prob = q_step2)
}
            