# created: 06/09/20
# modified: 06/09/20
# first trial of implementing simulation 1
# assume unit nonresponse follows ICIN, depends on weight
# and item nonresponse follow MNAR

V_adj = 1

# libraries
library(mvtnorm)
library(survey)

load("~/Missing/cps2012.rda")
dat12 <- da36383.0001
weights <- dat12$HWHHWGT/10000
# there are some zero weights entries (vacant)
nonzero_W <- weights[which(weights != 0)]

# Try equal weights, SRS
# nonzero_W <- rep(15,length(nonzero_W))
# nonzero_W <- nonzero_W*100

gen_pop_unit <- function(N=length(nonzero_W),HTW = nonzero_W,alpha = c(0.3,-0.5),
                         gamma = c(-0.25,0.1,-1.5),beta=c(0.5,-0.15),nu = c(-1.2,0.3)){
  # pi <- rbeta(N,2,5) 
  W <- HTW
  pi <- 1/W
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
  
  # sample U
  U <- X1
  pi_U <- pnorm(model.matrix(U~Y)%*%nu)
  U <- rbinom(N,1,p = pi_U)
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pi,U)))
}

getSample_Poisson <- function(population = pop_dat_Poisson){
  N <- dim(population)[1]
  pi <- population$pi*100  # increase inclusion probability 100 times to get enough sample
  getSampled <- rbinom(N,1,pi)
  sub_dat <- population[getSampled==1,]
  return(sub_dat)
}

# update Y and X per iteration
ANWC_unit_full2_1<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, 
                         gamma = gamma_s, nu = nu_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT,
                         pop_meanY = pop_meanY_HT, pop_sdY = pop_sdY_HT,weight_opt = "mean"){
  X1 <- sub_dat$X1
  W <- sub_dat$W
  Rx <- sub_dat$Rx
  Y <- sub_dat$Y
  U <- sub_dat$U
  
  # for U=1, make Rx unknown, so that we do not use them in Rx model
  Rx[which(U==1)] <- 999
  
  # for unit nonresponse, impute W based on weight option
  n_unit <- sum(U)
  n_mis <- sum(Rx==1)
  if(weight_opt == "mean"){
    W[which(U==1)] <- mean(W[which(U==0)])  
  }
  
  if(weight_opt == "max"){
    W[which(U==1)] <- max(W[which(U==0)])  
  }
  
  if(weight_opt == "bs"){
    W[which(U==1)] <- sample(W[which(U == 0)], n_unit, replace = TRUE)
  }
  
  # for U = 1, impute Y first
  pbeta <- 1/W
  m_Y <- glm(Y~pi,data = sub_dat[which(U == 0),],family=binomial(probit))
  p_Y <- predict(m_Y,newdata = data.frame(pi=pbeta[which(U == 1)]),type = "response")
  Y[which(U == 1)] <- rbinom(n_unit,1,prob = p_Y)
  
  # for U = 1, impute X
  m_X1 <- glm(X1~Y,data = sub_dat[which(U == 0),],family=binomial(probit))
  p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(U == 1)]),type = "response")
  X1[which(U == 1)] <- rbinom(n_unit,1,prob = p_X1)
  
  b0_alpha <- alpha; b0_alpha[] <- 0
  sigma0_alpha <- diag(alpha); diag(sigma0_alpha) <- 1
  b0_gamma <- gamma; b0_gamma[] <- 0
  sigma0_gamma <- diag(gamma); diag(sigma0_gamma) <- 1
  beta_Y <- t(rnorm(2))
  sigma0_beta <- diag(c(1,1))
  b0_beta <- rnorm(2)
  b0_nu <- nu; b0_nu[] <- 0
  sigma0_nu <- diag(nu); diag(sigma0_nu) <- 1
  
  alpha <- t(alpha_s)
  gamma <- t(gamma_s)
  nu <- t(nu_s)
  
  Y_accepted <- X1_accepted <- 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  BETA <- matrix(NA,niter-burnin,length(beta_Y))
  NU <- matrix(NA,niter-burnin,length(nu))
  X1_impute <- matrix(NA,niter-burnin,length(sub_dat$X1))
  Y_impute <- matrix(NA,niter-burnin,length(sub_dat$Y))
  # Rx_impute <- matrix(NA,niter-burnin,length(sub_dat$Rx))
  
  #GAMMA <- matrix(NA,niter,length(gamma))
  #ALPHA <- matrix(NA,niter,length(alpha))
  #BETA <- matrix(NA,niter,length(beta_Y))
  
  
  
  # prior on thetat
  theta0 <- c(1,1,1)
  
  for (i in 1:niter){
    # update alpha, with data augumentation
    # update Z|X,alpha
    Y_mat <- model.matrix(X1~as.factor(Y)+pbeta) # include weight in Y model
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
    YX_mat <- model.matrix(Rx~as.factor(Y)+X1,data = sub_dat[which(U==0),])
    G <- Rx[which(U==0)]
    G_mean <- YX_mat%*%t(gamma)
    G_mean_n <- G_mean[which(sub_dat[which(U==0),]$Rx == 0)] # Rx is based on all the data
    G_mean_p <- G_mean[which(sub_dat[which(U==0),]$Rx == 1)] # what we need is Rx==0 or 1for U==0 only
    U_Z <- G_mean; U_Z[] <- 0
    U_Z[which(sub_dat[which(U==0),]$Rx == 0)] <- runif(sum(Rx==0),pnorm(-Inf-G_mean_n),pnorm(0-G_mean_n))
    U_Z[which(sub_dat[which(U==0),]$Rx == 1)] <- runif(sum(Rx==1),pnorm(0-G_mean_p),pnorm(Inf-G_mean_p))
    G <- G_mean + qnorm(U_Z)
    
    # update gamma|G,Y,X
    sigma_hat <- solve(solve(sigma0_gamma) + t(YX_mat)%*%YX_mat)
    gamma_hat <- sigma_hat%*%(t(YX_mat)%*%G + solve(sigma0_gamma)%*%b0_gamma)
    gamma <- rmvnorm(1, mean = gamma_hat, sigma = sigma_hat)
    
    # for Rx = 1, and U == 0, sample X1_mis*
    Y_mat <- model.matrix(X1[which(U==0)]~as.factor(Y[which(U==0)])+pbeta[which(U==0)]) # use Y[which(U==0)] only 
    pr_X1_miss <- matrix(0,ncol=2,nrow=n_mis)
    colnames(pr_X1_miss) <- c("0","1")
    pi_x1 <- pnorm(0,Y_mat %*% t(alpha),1) 
    YX_matX0 <-  YX_matX1 <- YX_mat
    YX_matX0[,"X1"] <- rep(0,dim(YX_mat)[1])
    YX_matX1[,"X1"] <- rep(1,dim(YX_mat)[1])
    pi_Rx_temp0 <- dnorm(G,YX_matX0%*%t(gamma),1)
    pi_Rx_temp1<- dnorm(G,YX_matX1%*%t(gamma),1)
    pr_X1_miss[,"0"] <- (pi_Rx_temp0*pi_x1)[which(sub_dat[which(U==0),]$Rx == 1)]
    pr_X1_miss[,"1"] <- (pi_Rx_temp1*(1-pi_x1))[which(sub_dat[which(U==0),]$Rx == 1)]
    pr_X1_miss <- pr_X1_miss/matrix(rowSums(pr_X1_miss),ncol=2,nrow=n_mis)
    Ran_unif_X1_miss <- runif(nrow(pr_X1_miss))
    cumul_X1_miss <- pr_X1_miss%*%upper.tri(diag(ncol(pr_X1_miss)),diag=TRUE)
    
    proposed_X1 <- X1
    proposed_X1[which(sub_dat[which(U==0),]$Rx == 1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step for X1
    cur_mean <- sum(X1*W)
    proposed_mean <- sum(proposed_X1*W)
    log.r1 = dnorm(proposed_mean,mean = pop_mean,sd = pop_sd,log = TRUE) -  dnorm(cur_mean,mean = pop_mean,sd = pop_sd, log = TRUE)
    if(log(runif(1))< log.r1){
      X1 <- proposed_X1
      X1_accepted = X1_accepted + 1
    }
    
    # Update beta
    # update Q|beta,W
    # beta_Y <- coef(m_Y)
    Y_mat <- model.matrix(Y~pbeta) # include weight in Y model
    Q <- Y
    Q_mean <- Y_mat%*%t(beta_Y)
    Q_mean_n <- Q_mean[which(Y==0)]
    Q_mean_p <- Q_mean[which(Y==1)]
    U_Y <- Q_mean; U_Y[] <- 0
    U_Y[Y==0] <- runif(sum(Y==0),pnorm(-Inf-Q_mean_n),pnorm(0-Q_mean_n))
    U_Y[Y==1] <- runif(sum(Y==1),pnorm(0-Q_mean_p),pnorm(Inf-Q_mean_p))
    Q <- Q_mean + qnorm(U_Y)
    
    # update beta|Q,W
    sigma_hat <- solve(solve(sigma0_beta) + t(Y_mat)%*%Y_mat)
    beta_hat <- sigma_hat%*%(t(Y_mat)%*%Q + solve(sigma0_beta)%*%b0_beta)
    beta_Y <- rmvnorm(1, mean = beta_hat, sigma = sigma_hat)
    
    
    # Update nu
    # update V|nu,W
    U_mat <- model.matrix(U~Y) 
    V <- U
    V_mean <- U_mat%*%t(nu)
    V_mean_n <- V_mean[which(U==0)]
    V_mean_p <- V_mean[which(U==1)]
    V_U <- V_mean; V_U[] <- 0
    V_U[U==0] <- runif(sum(U==0),pnorm(-Inf-V_mean_n),pnorm(0-V_mean_n))
    V_U[V==1] <- runif(sum(U==1),pnorm(0-V_mean_p),pnorm(Inf-V_mean_p))
    V <- V_mean + qnorm(V_U)
    
    # update nu|V,W
    sigma_hat <- solve(solve(sigma0_nu) + t(U_mat)%*%U_mat)
    nu_hat <- sigma_hat%*%(t(U_mat)%*%V + solve(sigma0_nu)%*%b0_nu)
    nu <- rmvnorm(1, mean = nu_hat, sigma = sigma_hat)
    
    # for U = 1, sample Y_mis*
    pr_Y_miss <- matrix(0,ncol=2,nrow=n_unit)
    colnames(pr_Y_miss) <- c("0","1")
    pi_Y <- pnorm(0,Y_mat %*% t(beta_Y),1) 
    U_matY0 <-  U_matY1 <- U_mat
    U_matY0[,"Y"] <- rep(0,dim(U_mat)[1])
    U_matY1[,"Y"] <- rep(1,dim(U_mat)[1])
    pi_U_temp0 <- dnorm(U,U_matY0%*%t(nu),1)
    pi_U_temp1<- dnorm(U,U_matY1%*%t(nu),1)
    pr_Y_miss[,"0"] <- (pi_U_temp0*pi_Y)[which(U==1)]
    pr_Y_miss[,"1"] <- (pi_U_temp1*(1-pi_Y))[which(U==1)]
    pr_Y_miss <- pr_Y_miss/matrix(rowSums(pr_Y_miss),ncol=2,nrow=n_unit)
    Ran_unif_Y_miss <- runif(nrow(pr_Y_miss))
    cumul_Y_miss <- pr_Y_miss%*%upper.tri(diag(ncol(pr_Y_miss)),diag=TRUE)
    
    # Update Y
    proposed_Y <- Y
    proposed_Y[which(U==1)] <- rowSums(Ran_unif_Y_miss > cumul_Y_miss)
    
    # M-H step for Y
    cur_meanY <- sum(Y*W) 
    proposed_meanY <- sum(proposed_Y*W)
    log.r2 = dnorm(proposed_meanY,mean = pop_meanY,sd = pop_sdY,log = TRUE) -  dnorm(cur_meanY,mean = pop_meanY,sd = pop_sdY, log = TRUE)
    if(log(runif(1))< log.r2){
      Y <- proposed_Y
      Y_accepted = Y_accepted + 1
    }
    
    
    # Update X where U == 1
    p_X1 <- pnorm(model.matrix(X1~as.factor(Y)+pbeta)%*%t(alpha))
    X1[which(U == 1)] <- rbinom(n_unit,1,prob = p_X1[which(U == 1)])
    
    
    # resample missing weight if using bootstrap
    if(weight_opt == "bs"){
      W[which(U==1)] <- sample(W[which(U == 0)], n_unit, replace = TRUE)
    }
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      BETA[i-burnin,] <- beta_Y
      NU[i-burnin,] <- nu
      X1_impute[i-burnin,] <- X1
      Y_impute[i-burnin,] <- Y
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  YL <- Y_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA, beta = BETA,nu = NU,MI_dataX = XL, MI_dataY = YL, acceptRatio = X1_accepted/niter))
}

getResults_unit2_1 <- function(dataMI_X,dataMI_Y,alpha_len, gamma_len,beta_len, nu_len, sub_dat){
  n <- dim(dataMI_X)[1]
  ans <- matrix(NA,n,alpha_len)
  ul <- matrix(NA,n,alpha_len)
  total_var <- total_varY <- rep(NA,n)
  
  # calculate beta
  ans_b <- matrix(NA,n,beta_len)
  ul_b <- matrix(NA,n,beta_len)
  
  # calculate nu
  ans_nu <- matrix(NA,n,nu_len)
  ul_nu <- matrix(NA,n,nu_len)
  
  # calculate gamma using stats package
  ans_g <- matrix(NA,n,gamma_len)
  ul_g <- matrix(NA,n,gamma_len)
  
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI_X[i,],rownames(sub_dat),sub_dat$W,dataMI_Y[i,],sub_dat$Rx,sub_dat$U),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y","Rx","U")
    total_var[i] <- sum((as.numeric(test_dat$X1)*sub_dat$W)^2*(1-1/sub_dat$W))
    total_varY[i] <- sum((as.numeric(test_dat$Y)*sub_dat$W)^2*(1-1/sub_dat$W))
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$pi <- 1/test_dat$W
    test_dat$Y <- as.factor(test_dat$Y)
    test_dat$Rx <- as.factor(test_dat$Rx)
    test_dat$U <- as.factor(test_dat$U)
    test_dat2 <- test_dat[which(test_dat$U == 0),]
    mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
    mydesign2 <- svydesign(id = ~id,data = test_dat2,weight = ~W)
    m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
    ans[i,] <- coef(m1)
    ul[i,] <- diag(vcov(m1))
    m3 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
    m2 <- svyglm(Rx~Y+X1,design = mydesign2,family = quasibinomial(link = "probit"))
    ans_b[i,] <- coef(m3)
    ul_b[i,] <- diag(vcov(m3))
    ans_g[i,] <- coef(m2)
    ul_g[i,] <- diag(vcov(m2))
    m4 <- svyglm(U~Y,design = mydesign,family = quasibinomial(link = "probit"))
    ans_nu[i,] <- coef(m4)
    ul_nu[i,] <- diag(vcov(m4))
  }
  # colMeans(ans alpha)
  u_L_bar <- colMeans(ul)
  b_L <- apply((scale(ans,scale=FALSE))^2/(n-1),2,sum)
  T_L <- (1+1/n)*b_L + u_L_bar
  
  # colMeans(ans_g)
  u_L_g <- colMeans(ul_g)
  b_L_g <- apply((scale(ans_g,scale=FALSE))^2/(n-1),2,sum)
  T_L_g <- (1+1/n)*b_L_g + u_L_g
  
  # colMeans(ans_b)
  u_L_b <- colMeans(ul_b)
  b_L_b <- apply((scale(ans_b,scale=FALSE))^2/(n-1),2,sum)
  T_L_b <- (1+1/n)*b_L_b + u_L_b
  
  # colMeans(ans_nu)
  u_L_nu <- colMeans(ul_nu)
  b_L_nu <- apply((scale(ans_nu,scale=FALSE))^2/(n-1),2,sum)
  T_L_nu <- (1+1/n)*b_L_nu + u_L_nu
  
  total <- total_Y <- rep(NA,n)
  for (i in 1:n){
    total[i] <- sum(dataMI_X[i,]*sub_dat$W)
    total_Y[i] <- sum(dataMI_Y[i,]*sub_dat$W)
  }
  u_L_t <- mean(total_var)
  b_L_t <- apply((scale(total,scale=FALSE))^2/(n-1),2,sum)
  T_L_t <- (1+1/n)*b_L_t + u_L_t
  u_L_tY <- mean(total_varY)
  b_L_tY <- apply((scale(total_Y,scale=FALSE))^2/(n-1),2,sum)
  T_L_tY <- (1+1/n)*b_L_tY + u_L_tY
  return(list(alpha = ans, gamma = ans_g, beta = ans_b, nu = ans_nu, total_meanX = mean(total), total_varX = T_L_t,
              total_meanY = mean(total_Y), total_varY = T_L_tY,
              alpha_var = T_L, alpha_mean = colMeans(ans),gamma_mean = colMeans(ans_g),beta_mean = colMeans(ans_b),
              nu_mean = colMeans(ans_nu), gamma_var = T_L_g,beta_var = T_L_b, nu_var = T_L_nu))
}


sub_dat <- getSample_Poisson(population = pop_dat)
sub_dat$W <- sub_dat$W/100
sub_dat$pi <- sub_dat$pi*100

# HT estimator of population total
X1 <- sub_dat$X1
W <- sub_dat$W
Rx <- sub_dat$Rx
Y <- sub_dat$Y
U <- sub_dat$U
pi <- sub_dat$pi
pbeta <- 1/sub_dat$W
n_mis <- sum(sub_dat$Rx == 1)
m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0 & U == 0),],family=binomial(probit))
p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1||U == 1)]),type = "response")
X1[which(Rx == 1 || U == 1)] <- rbinom(length(p_X1),1,prob = p_X1)
pop_sd_HT <- sqrt(sum((X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj 
pop_mean_HT <- sum(pop_dat$X1)

m_Y <- glm(Y~pi,data = sub_dat[which(U == 0),],family=binomial(probit))
p_Y <- predict(m_Y,newdata = data.frame(pi=pi[which(U == 1)]),type = "response")
Y[which(U == 1)] <- rbinom(length(p_Y),1,prob = p_Y)
pop_sdY_HT <- sqrt(sum((Y/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj 
pop_meanY_HT <- sum(pop_dat$Y)

alpha_s <- rnorm(3)
gamma_s <- rnorm(3)
nu_s <- rnorm(2)
testListANWC <- ANWC_unit_full2_1(weight_opt = "known")
resultListANWC <- getResults_unit2_1(dataMI_X = testListANWC$MI_dataX, dataMI_Y = testListANWC$MI_dataY,
                                  alpha_len=2, gamma_len = 3,beta_len = 2, nu_len = 2, sub_dat)

ratio <- testListANWC$acceptRatio
MI_dataANWC <- data.frame(X = testListANWC$MI_dataX, Y = testListANWC$MI_dataY)

save(resultListANWC,sub_dat,ratio,MI_dataANWC,pop_sd_HT,file = paste("./unit6/Mis_",run_index,".RData",sep=""))
