# created: 05/08
# modified: 05/12
# unit nonresponse final trial to be run on server
# update Y and X per iteration, too
# do not update Rx per iteration
# do not impute missing X for U==1
# fully Bayesian
# User can choose how to deal with unit nonresponse weights

load("./unit_data1.RData")
# load("./unit_SRS.RData")

run_index = commandArgs(trailingOnly = TRUE)
V_adj = 1

# libraries
library(mvtnorm)
library(survey)

getSample_Poisson <- function(population = pop_dat_Poisson){
  N <- dim(population)[1]
  pi <- population$pi*100  # increase inclusion probability 100 times to get enough sample
  getSampled <- rbinom(N,1,pi)
  sub_dat <- population[getSampled==1,]
  return(sub_dat)
}

# update Y and X per iteration
ANWC_unit_full<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, 
                         gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT, weight_opt = "mean"){
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
  
  alpha <- t(alpha_s)
  gamma <- t(gamma_s)
  
  X1_accepted = 0
  
  GAMMA <- matrix(NA,niter-burnin,length(gamma))
  ALPHA <- matrix(NA,niter-burnin,length(alpha))
  BETA <- matrix(NA,niter-burnin,length(beta_Y))
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
    proposed_X1[which(Rx==1)] <- rowSums(Ran_unif_X1_miss>cumul_X1_miss)
    
    # M-H step
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
    
    # Update Y
    p_Y <- pnorm(beta_Y[1]+beta_Y[2]*pbeta)
    Y[which(U == 1)] <- rbinom(n_unit,1,prob = p_Y[which(U == 1)])
    
    
    # Update X
    p_X1 <- pnorm(model.matrix(X1~as.factor(Y)+pbeta)%*%t(alpha))
    X1[which(U == 1)] <- rbinom(n_unit,1,prob = p_X1[which(U == 1)])
    
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      BETA[i-burnin,] <- beta_Y
      X1_impute[i-burnin,] <- X1
      Y_impute[i-burnin,] <- Y
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  YL <- Y_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA, beta = BETA,MI_dataX = XL, MI_dataY = YL, acceptRatio = X1_accepted/niter))
}

getResults_unit <- function(dataMI_X,dataMI_Y,alpha_len, gamma_len,beta_len, sub_dat){
  n <- dim(dataMI_X)[1]
  ans <- matrix(NA,n,alpha_len)
  ul <- matrix(NA,n,alpha_len)
  total_var <- rep(NA,n)
  
  # calculate beta
  ans_b <- matrix(NA,n,beta_len)
  ul_b <- matrix(NA,n,beta_len)
  
  # calculate gamma using stats package
  ans_g <- matrix(NA,n,gamma_len)
  ul_g <- matrix(NA,n,gamma_len)
  
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI_X[i,],rownames(sub_dat),sub_dat$W,dataMI_Y[i,],sub_dat$Rx,sub_dat$U),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y","Rx","U")
    total_var[i] <- sum((as.numeric(test_dat$X1)*sub_dat$W)^2*(1-1/sub_dat$W))
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
    
  }
  # colMeans(ans)
  u_L_bar <- colMeans(ul)
  b_L <- apply((scale(ans,scale=FALSE))^2/(n-1),2,sum)
  T_L <- (1+1/n)*b_L + u_L_bar
  
  # colMeans(ans_g)
  u_L_g <- colMeans(ul_g)
  b_L_g <- apply((scale(ans_g,scale=FALSE))^2/(n-1),2,sum)
  T_L_g <- (1+1/n)*b_L_g + u_L_g
  
  # colMeans(ans_g)
  u_L_b <- colMeans(ul_b)
  b_L_b <- apply((scale(ans_b,scale=FALSE))^2/(n-1),2,sum)
  T_L_b <- (1+1/n)*b_L_b + u_L_b
  
  total <- rep(NA,n)
  for (i in 1:n){
    total[i] <- sum(dataMI_X[i,]*sub_dat$W)
  }
  u_L_t <- mean(total_var)
  b_L_t <- apply((scale(total,scale=FALSE))^2/(n-1),2,sum)
  T_L_t <- (1+1/n)*b_L_t + u_L_t
  return(list(alpha = ans, gamma = ans_g, beta = ans_b, total_mean = mean(total), total_var = T_L_t, 
              alpha_var = T_L, alpha_mean = colMeans(ans),gamma_mean = colMeans(ans_g),beta_mean = colMeans(ans_b),
              gamma_var = T_L_g,beta_var = T_L_b))
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
pbeta <- 1/sub_dat$W
n_mis <- sum(sub_dat$Rx == 1)
m_X1 <- glm(X1~Y,data = sub_dat[which(Rx ==0 & U == 0),],family=binomial(probit))
p_X1 <- predict(m_X1,newdata = data.frame(Y=Y[which(Rx == 1||U == 1)]),type = "response")
X1[which(Rx == 1 || U == 1)] <- rbinom(length(p_X1),1,prob = p_X1)
pop_sd_HT <- sqrt(sum((X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj 
pop_mean_HT <- sum(pop_dat$X1)

alpha_s <- rnorm(3)
gamma_s <- rnorm(3)
testListANWC <- ANWC_unit_full(weight_opt = "known")
resultListANWC <- getResults_unit(dataMI_X = testListANWC$MI_dataX, dataMI_Y = testListANWC$MI_dataY,
                                  alpha_len=2, gamma_len = 3,beta_len = 2, sub_dat)

ratio <- testListANWC$acceptRatio
MI_dataANWC <- data.frame(X = testListANWC$MI_dataX, Y = testListANWC$MI_dataY)

save(resultListANWC,sub_dat,ratio,MI_dataANWC,pop_sd_HT,file = paste("./unit9/Mis_",run_index,".RData",sep=""))
