gen_pop_Poisson <- function(N=50000,N1=35000, N2= 15000,alpha = c(0.3,-0.5),gamma = c(-0.25,0.1,-1.1),beta=c(0.5,-0.15)){
  pi <- rbeta(N,2,5) 
  W <- 1/pi
  Y <- W
  
  # generate Y give weights, have some pi_Y very close to 0
  pi_Y <- pnorm(model.matrix(Y~pi)%*%beta)
  Y <- rbinom(N,1,p=pi_Y)
  
  # generate X1 given Y
  X1 <- Y
  pi_x1 <- pnorm(model.matrix(X1~Y)%*%alpha) # pnorm is CDF or normal
  X1 <- rbinom(N,1,p=pi_x1)
  
  # sample Rx|X,Y
  Rx <- X1
  pi_Rx <- pnorm(model.matrix(Rx~Y+X1)%*%gamma)
  Rx <- rbinom(N,1,p = pi_Rx)
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pi)))
}

gen_pop_PPS <- function(N=50000,alpha = c(0.5,-0.5,-1),gamma = c(-0.6,0.1,0.3,-1.1),beta= c(2,-1,0.8),theta = c(0.1,0.45,0.45)){
  alpha0 <- alpha[1]
  alpha12 <- alpha[2]
  alpha13 <- alpha[3]
  gamma0 <- gamma[1]
  gamma12 <- gamma[2]
  gamma13 <- gamma[3]
  gamma2 <- gamma[4]
  beta12 <- beta[1]
  beta13 <- beta[2]
  beta2 <- beta[3]
  
  # generate Y
  Y <- sample(c(1,2,3),N,replace = TRUE, prob = theta)
  # Y[1:N1] <- sample(c(1,2,3),N1,replace = TRUE, prob = theta1)
  # Y[(N1+1):N] <- sample(c(1,2,3),N2,replace = TRUE, prob = theta2)
  
  # generate X1 given Y
  pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) # pnorm is CDF or normal
  X1 <- rbinom(N,1,p=pi_x1)
  
  # sample Rx|X,Y
  pi_Rx <- pnorm(gamma0 + gamma12*(Y == 2) + gamma13*(Y==3) + gamma2*X1)
  Rx <- rbinom(N,1,p = pi_Rx)
  
  # generate Z|X,Y, use this variable for weight
  Z <- rnorm(1,-2.5,1) + beta12*(Y==2) + beta13*(Y==3) + beta2*X1
  
  # generate inclusion probability
  pi <- pnorm(Z)
  W <- 1/pi
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pi)))
}

gen_pop <- function(N=50000,N1=35000, N2= 15000, n1 = 1500, n2 = 3500, alpha = c(0.5,-0.5,-1),gamma = c(-0.25,0.1,0.3,-1.1),theta1= c(0.5,0.15,0.35),theta2 = c(0.1,0.45,0.45)){
  alpha0 <- alpha[1]
  alpha12 <- alpha[2]
  alpha13 <- alpha[3]
  gamma0 <- gamma[1]
  gamma12 <- gamma[2]
  gamma13 <- gamma[3]
  gamma2 <- gamma[4]
  Y <- rep(NA,N)
  Y[1:N1] <- sample(c(1,2,3),N1,replace = TRUE, prob = theta1)
  Y[(N1+1):N] <- sample(c(1,2,3),N2,replace = TRUE, prob = theta2)
  W <- c(rep(N1/n1,N1),rep(N2/n2,N2))
  pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) # pnorm is CDF or normal
  
  # generate X1 given Y
  X1 <- rbinom(N,1,p=pi_x1)
  
  # sample Rx|X,Y
  pi_Rx <- pnorm(gamma0 + gamma12*(Y == 2) + gamma13*(Y==3) + gamma2*X1)
  Rx <- rbinom(N,1,p = pi_Rx)
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W)))
}
# call it
getSample_Poisson <- function(population = pop_dat_Poisson,n=5000,N=50000){
  getSampled <- rbinom(N,1,population$pi)
  sub_dat <- population[getSampled==1,]
  return(sub_dat)
}

getSampleStatsPoisson <- function(population = pop_dat_Poisson,N=50000){
  getSampled <- rbinom(N,1,population$pi)
  sub_dat <- population[getSampled==1,]
  return(sum(sub_dat$X1*sub_dat$W))
}

getSample <- function(population,n1=1500,n2=3500,N1=35000, N2= 15000){
  getSampled <- c(sample(1:N1,n1),sample((N1+1):(N1+N2),n2))
  sub_dat <- population[getSampled,]
  return(sub_dat)
}

getSampleStats <- function(population,n1=1500,n2=3500,N1=35000, N2= 15000){
  getSampled <- c(sample(1:N1,n1),sample((N1+1):(N1+N2),n2))
  sub_dat <- population[getSampled,]
  return(sum(sub_dat$X1*sub_dat$W))
}

doGibbs4<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  sub_dat$X1 <- as.factor(sub_dat$X1)
  sub_dat$Y <- as.factor(sub_dat$Y)
  
  Model <- list(X1=~Y,Rx=~Y+X1)
  alpha <- rnorm(ncol(model.matrix(Model$X1,data=sub_dat)))
  gamma <- rnorm(ncol(model.matrix(Model$R,data=sub_dat)))
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  Z <- data.frame(X1=rnorm(nrow(sub_dat)))
  Z$R <- rnorm(nrow(sub_dat))
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    
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
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_Rx_temp0 <- dnorm(Z$R,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    pi_Rx_temp1 <- dnorm(Z$R,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbs<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # for Rx = 1, sample X1_mis*
    # pi_x1 <- pnorm(alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3)) 
    # pi_Rx_temp <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3)
    # q1 <- pi_x1/(1-pi_x1)*pnorm(pi_Rx_temp+gamma[4])/pnorm(pi_Rx_temp)
    # proposed_X1 <- X1
    # proposed_X1[which(Rx==1)] <- rbinom(n_mis, 1, p = q1/(1+q1))
    
    # M-H step
    # cur_mean <- sum(X1*W)
    # proposed_mean <- sum(proposed_X1*W)
    # log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log=TRUE) - dnorm(cur_mean,mean = pop_mean,sd = pop_sd,log=TRUE)
    # if(log(runif(1))< log.r1){
    #  X1 <- proposed_X1
    #  X1_accepted = X1_accepted + 1
    # }
    
    # update alpha, with data augumentation
    # update Z|X,alpha
    Z <- X1
    Z_mean <- alpha[1] + alpha[2]*(Y==2) + alpha[3]*(Y==3)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    # Y_n <- Y[which(X1 == 0)]
    # Z_mean_n <- alpha[1] + alpha[2]*(Y_n==2) + alpha[3]*(Y_n==3)
    # Z[which(X1 == 0)] <- rtruncnorm(1, a=-Inf, b=0, mean = Z_mean_n, sd = 1)
    # Y_p <- Y[which(X1 != 0)]
    # Z_mean_p <- alpha[1] + alpha[2]*(Y_p==2) + alpha[3]*(Y_p==3)
    #Z[which(X1 != 0)] <- rtruncnorm(1, a=0, b=Inf, mean = Z_mean_p, sd = 1)
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    Y_mat <- model.matrix(X1~as.factor(Y))
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    G <- Rx
    G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4]*X1
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    #Y_n <- Y[which(Rx == 0)]
    #X_n <- X1[which(Rx == 0)]
    #G_mean_n <- gamma[1] + gamma[2]*(Y_n==2) + gamma[3]*(Y_n==3) + gamma[4]*X_n
    #G[which(Rx == 0)] <- rtruncnorm(1, a=-Inf, b=0, mean = G_mean_n, sd = 1)
    #Y_p <- Y[which(Rx != 0)]
    #X_p <- X1[which(Rx != 0)]
    #G_mean_p <- gamma[1] + gamma[2]*(Y_p==2) + gamma[3]*(Y_p=3) + gamma[4]*X_p
    #G[which(Rx != 0)] <- rtruncnorm(1, a=0, b=Inf, mean = G_mean_p, sd = 1)
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsGoodStart<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Z <- X1
    Z_mean <- alpha[1] + alpha[2]*(Y==2) + alpha[3]*(Y==3)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    # Y_n <- Y[which(X1 == 0)]
    # Z_mean_n <- alpha[1] + alpha[2]*(Y_n==2) + alpha[3]*(Y_n==3)
    # Z[which(X1 == 0)] <- rtruncnorm(1, a=-Inf, b=0, mean = Z_mean_n, sd = 1)
    # Y_p <- Y[which(X1 != 0)]
    # Z_mean_p <- alpha[1] + alpha[2]*(Y_p==2) + alpha[3]*(Y_p==3)
    #Z[which(X1 != 0)] <- rtruncnorm(1, a=0, b=Inf, mean = Z_mean_p, sd = 1)
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    Y_mat <- model.matrix(X1~as.factor(Y))
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    G <- Rx
    G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4]*X1
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    #Y_n <- Y[which(Rx == 0)]
    #X_n <- X1[which(Rx == 0)]
    #G_mean_n <- gamma[1] + gamma[2]*(Y_n==2) + gamma[3]*(Y_n==3) + gamma[4]*X_n
    #G[which(Rx == 0)] <- rtruncnorm(1, a=-Inf, b=0, mean = G_mean_n, sd = 1)
    #Y_p <- Y[which(Rx != 0)]
    #X_p <- X1[which(Rx != 0)]
    #G_mean_p <- gamma[1] + gamma[2]*(Y_p==2) + gamma[3]*(Y_p=3) + gamma[4]*X_p
    #G[which(Rx != 0)] <- rtruncnorm(1, a=0, b=Inf, mean = G_mean_p, sd = 1)
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsSwitch<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Z <- X1
    Z_mean <- alpha[1] + alpha[2]*(Y==2) + alpha[3]*(Y==3)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    Y_mat <- model.matrix(X1~as.factor(Y))
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    G <- Rx
    G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4]*X1
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    proposed_X1 <- X1
    proposed_X1[sample(which(Rx == 1),1)] <- 1 - proposed_X1[sample(which(Rx == 1),1)]
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsSwitchGoodStart<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  # X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Z <- X1
    Z_mean <- alpha[1] + alpha[2]*(Y==2) + alpha[3]*(Y==3)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    Y_mat <- model.matrix(X1~as.factor(Y))
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    G <- Rx
    G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4]*X1
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    proposed_X1 <- X1
    proposed_X1[sample(which(Rx == 1),1)] <- 1 - proposed_X1[sample(which(Rx == 1),1)]
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsMarWeight<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  # alpha should be length 4 
  # gamma can be unchanged
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  S <- ifelse(sub_dat$W >5,1,2)
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Z <- X1
    Z_mean <- alpha[1] + alpha[2]*(Y==2) + alpha[3]*(Y==3) + alpha[4]*(S==2)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    Y_mat <- model.matrix(X1~as.factor(Y)+S)
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    G <- Rx
    G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3)
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
    # cur_mean <- sum(X1*W)
    # proposed_mean <- sum(proposed_X1*W)
    # log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    # if(log(runif(1))< log.r1){
    #  X1 <- proposed_X1
    #  X1_accepted = X1_accepted + 1
    # }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsMarWeightPoisson<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  pbeta <- 1/sub_dat$W
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  alpha <- t(alpha_s)
  gamma <- t(gamma_s)
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Y_mat <- model.matrix(X1~as.factor(Y)+pbeta)
    Z <- X1
    #Z_mean <- alpha[1] + alpha[2]*(Y==2) + alpha[3]*(Y==3) + alpha[4]*(1/W)
    Z_mean <- Y_mat%*%t(alpha)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    G <- Rx
    # G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3)
    G_mean <- YX_mat%*%t(gamma)
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    # pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_x1 <- pnorm(0,Y_mat %*% t(alpha),1) 
    YX_matX0 <-  YX_matX1 <- YX_mat
    YX_matX0[,"X1"] <- rep(0,dim(YX_mat)[1])
    YX_matX1[,"X1"] <- rep(1,dim(YX_mat)[1])
    #pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    #pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pi_Rx_temp0 <- dnorm(G,YX_matX0%*%t(gamma),1)
    pi_Rx_temp1<- dnorm(G,YX_matX1%*%t(gamma),1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
    # cur_mean <- sum(X1*W)
    # proposed_mean <- sum(proposed_X1*W)
    # log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    # if(log(runif(1))< log.r1){
    #  X1 <- proposed_X1
    #  X1_accepted = X1_accepted + 1
    # }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsANWCPoisson<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  pbeta <- 1/sub_dat$W
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  alpha <- t(alpha_s)
  gamma <- t(gamma_s)
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Y_mat <- model.matrix(X1~as.factor(Y)+pbeta)
    Z <- X1
    Z_mean <- Y_mat%*%t(alpha)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    G <- Rx
    # G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3)
    G_mean <- YX_mat%*%t(gamma)
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    # pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_x1 <- pnorm(0,Y_mat %*% t(alpha),1) 
    YX_matX0 <-  YX_matX1 <- YX_mat
    YX_matX0[,"X1"] <- rep(0,dim(YX_mat)[1])
    YX_matX1[,"X1"] <- rep(1,dim(YX_mat)[1])
    #pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    #pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pi_Rx_temp0 <- dnorm(G,YX_matX0%*%t(gamma),1)
    pi_Rx_temp1<- dnorm(G,YX_matX1%*%t(gamma),1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
    X1 <- proposed_X1
    X1_accepted = X1_accepted + 1
  }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

doGibbsSwitchWCPoisson<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  pbeta <- 1/sub_dat$W
  n_mis <- sum(sub_dat$Rx == 1)
  m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1)]),type = "response")
  X1[which(Rx == 1)] <- rbinom(n_mis,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  
  alpha <- t(alpha_s)
  gamma <- t(gamma_s)
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Y_mat <- model.matrix(X1~as.factor(Y)+pbeta)
    Z <- X1
    Z_mean <- Y_mat%*%t(alpha)
    Z_mean_n <- Z_mean[which(X1==0)]
    Z_mean_p <- Z_mean[which(X1==1)]
    U_X1 <- Z_mean; U_X1[] <- 0
    U_X1[X1==0] <- runif(sum(X1==0),pnorm(-Inf-Z_mean_n),pnorm(0-Z_mean_n))
    U_X1[X1==1] <- runif(sum(X1==1),pnorm(0-Z_mean_p),pnorm(Inf-Z_mean_p))
    Z <- Z_mean + qnorm(U_X1)
    
    # update alpha|Z,Y
    sigma_hat <- solve(solve(sigma0_alpha) + t(Y_mat)%*%Y_mat)
    alpha_hat <- sigma_hat%*%(t(Y_mat)%*%Z + solve(sigma0_alpha)%*%b0_alpha)
    alpha <- rmvnorm(1, mean = alpha_hat, sigma = sigma_hat)
    
    # update gamma
    # update G|Rx, gamma
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1)
    G <- Rx
    # G_mean <- gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3)
    G_mean <- YX_mat%*%t(gamma)
    G_mean_n <- G_mean[which(Rx==0)]
    G_mean_p <- G_mean[which(Rx==1)]
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[Rx==0] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[Rx==1] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update alpha|Z,Y
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, sample X1_mis*
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    # pi_x1 <- pnorm(0,alpha[1] + alpha[2]*(Y == 2) + alpha[3]*(Y==3),1) 
    pi_x1 <- pnorm(0,Y_mat %*% t(alpha),1) 
    YX_matX0 <-  YX_matX1 <- YX_mat
    YX_matX0[,"X1"] <- rep(0,dim(YX_mat)[1])
    YX_matX1[,"X1"] <- rep(1,dim(YX_mat)[1])
    #pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3),1)
    #pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*(Y==2) + gamma[3]*(Y==3) + gamma[4],1)
    pi_Rx_temp0 <- dnorm(G,YX_matX0%*%t(gamma),1)
    pi_Rx_temp1<- dnorm(G,YX_matX1%*%t(gamma),1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(Rx==1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(Rx==1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    # for Rx = 1, sample X1_mis*
    proposed_X1 <- X1
    proposed_X1[sample(which(Rx == 1),1)] <- 1 - proposed_X1[sample(which(Rx == 1),1)]
    
    # M-H step
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}
getResults <- function(dataMI,alpha_len, gamma_len,sub_dat){
  n <- dim(dataMI)[1]
  ans <- matrix(NA,n,alpha_len)
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI[i,],rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y")
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$Y <- as.factor(test_dat$Y)
    mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
    m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
    ans[i,] <- coef(m1)
  }
  colMeans(ans)
  # 0.5,-0.5,-1
  
  # calculate gamma using stats package
  ans_g <- matrix(NA,n,gamma_len)
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI[i,],rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$Rx),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y","Rx")
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$Y <- as.factor(test_dat$Y)
    test_dat$Rx <- as.factor(test_dat$Rx)
    mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
    m2 <- glm(Rx~Y+X1,data = test_dat, family = binomial(link = "probit"))
    ans_g[i,] <- coef(m2)
  }
  colMeans(ans_g)
  total <- rep(NA,n)
  for (i in 1:n){
    total[i] <- sum(dataMI[i,]*sub_dat$W)
  }
  return(list(alpha = ans, gamma = ans_g, pop_total =  total))
}

getResultsStrata <- function(dataMI, n1 = 1500, n2 = 3500,sub_dat){
  n <- dim(dataMI)[1]
  ans <- matrix(NA,n,3)
  for (i in 1:n){
    str <- c(rep(1,n1),rep(2,n2))
    test_dat <- as.data.frame(cbind(dataMI[i,],rownames(sub_dat),sub_dat$W,str,sub_dat$Y),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","str","Y")
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$str <- as.numeric(test_dat$str)
    test_dat$Y <- as.factor(test_dat$Y)
    mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W,strata = ~str)
    m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
    ans[i,] <- coef(m1)
  }
  colMeans(ans)
  # 0.5,-0.5,-1
  
  # calculate gamma using stats package
  ans_g <- matrix(NA,n,4)
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI[i,],rownames(sub_dat),sub_dat$W,str,sub_dat$Y,sub_dat$Rx),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","str","Y","Rx")
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
  total <- rep(NA,n)
  for (i in 1:n){
    total[i] <- sum(dataMI[i,]*sub_dat$W)
  }
  return(list(alpha = ans, gamma = ans_g, pop_total =  total))
}