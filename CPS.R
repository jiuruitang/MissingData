# use CPS2012 data weights to generate X and Y
# created: 01/20
# modified: 
load("~/Missing/cps2012.rda")
weights <- dat12$HWHHWGT/10000
# there are some zero weights entries (vacant)
nonzero_W <- weights[which(weights != 0)]

gen_pop_Poisson <- function(N=length(nonzero_W),HTW = nonzero_W,alpha = c(0.3,-0.5),gamma = c(-0.25,0.1,-1.5),beta=c(0.5,-0.15)){
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
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pi)))
}

getSample_Poisson <- function(population = pop_dat_Poisson){
  N <- dim(population)[1]
  pi <- population$pi*100  # increase inclusion probability 100 times to get enough sample
  getSampled <- rbinom(N,1,pi)
  sub_dat <- population[getSampled==1,]
  return(sub_dat)
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

#######################
### use functions
library(mvtnorm)
library(survey)
pop_dat <- gen_pop_Poisson() # 24.3% missing
sub_dat <- getSample_Poisson(population = pop_dat)
# adjust weight and pi accordingly
sub_dat$W <- sub_dat$W/100
sub_dat$pi <- sub_dat$pi*100
pop_mean_HT <- sum(pop_dat$X1) # population total truth: 
sum(sub_dat$X1*sub_dat$W) # best estimation based on full sample: 62758.76
# pop_sd_HT <- sqrt(sum((sub_dat$X1)^2/sub_dat$pi*(1-sub_dat$pi))) # HT estimator of population total
pop_sd_HT <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi))) 
# pop_sd_HT <- sqrt(sum(((pop_dat$X1/pop_dat$pi)^2/(pop_dat$pi))*(1-pop_dat$pi)))

alpha_s <- rnorm(3)
gamma_s <- rnorm(3)
testListANWC <- doGibbsANWCPoisson()
resultListANWC <- getResults(dataMI = testListANWC$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
colMeans(resultListANWC$alpha)
colMeans(resultListANWC$gamma)
mean(resultListANWC$pop_total)


### repeat 10 times
alpha_MAR <- matrix(NA,10,2)
gamma_MAR <- matrix(NA,10,3)
total_MAR <- rep(NA,10)
alpha_ANWC <- matrix(NA,10,2)
gamma_ANWC <- matrix(NA,10,3)
total_ANWC <- rep(NA,10)
acceptRatio_ANWC <- rep(NA,10)
total_truth <- rep(NA,10)
#alpha_R <- matrix(NA,10,2)
#gamma_R <- matrix(NA,10,3)
#total_R <- rep(NA,10)
#acceptRatio_R <- rep(NA,10)
alpha_truth <- matrix(NA,10,2)
gamma_truth <- matrix(NA,10,3)
for(i in 1:10){
  sub_dat <- getSample_Poisson(population = pop_dat)
  sub_dat$W <- sub_dat$W/100
  sub_dat$pi <- sub_dat$pi*100
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y+ I(1/W),design = mydesign,family = quasibinomial(link = "probit"))
  alpha_truth[i,] <- coef(m1)[1:2]
  m2 <- glm(Rx~Y+X1+I(1/W),data = sub_dat, family = binomial(link = "probit"))
  gamma_truth[i,] <- coef(m2)[1:3]
  
  alpha_s <- rnorm(3)
  gamma_s <- rnorm(3)
  testListMAR <- doGibbsMarWeightPoisson()
  resultListMAR <- getResults(dataMI = testListMAR$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
  alpha_MAR[i,] <- colMeans(resultListMAR$alpha)
  gamma_MAR[i,] <- colMeans(resultListMAR$gamma)
  total_MAR[i] <- mean(resultListMAR$pop_total)
  
  testListANWC <- doGibbsANWCPoisson()
  resultListANWC <- getResults(dataMI = testListANWC$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
  alpha_ANWC[i,] <- colMeans(resultListANWC$alpha)
  gamma_ANWC[i,] <- colMeans(resultListANWC$gamma)
  total_ANWC[i] <- mean(resultListANWC$pop_total)
  acceptRatio_ANWC[i] <- testListANWC$acceptRatio
  total_truth[i] <- sum(sub_dat$W*sub_dat$X1)
}

resultMat <- matrix(NA,6,4)
resultMat[,1] <-  c(colMeans(alpha_MAR),colMeans(gamma_MAR),mean(total_MAR))
resultMat[,2] <-  c(apply(alpha_MAR,2,sd)/sqrt(10),apply(gamma_MAR,2,sd)/sqrt(10),sd(total_MAR))
resultMat[,3] <-  c(colMeans(alpha_ANWC),colMeans(gamma_ANWC),mean(total_ANWC))
resultMat[,4] <-  c(apply(alpha_ANWC,2,sd)/sqrt(10),apply(gamma_ANWC,2,sd)/sqrt(10),sd(total_ANWC))
poisson_result <- as.data.frame(resultMat)
names(poisson_result) <- c("MAR+W Mean","MAR+W SD","AN+WC Mean", "AN+WC SD")

poisson_result$Sampletruth <- c(colMeans(alpha_truth),colMeans(gamma_truth),mean(total_truth))
poisson_result$truth <-c(c(0.3,-0.5),c(-0.25,0.1,-1.5),sum(pop_dat$X1))
poisson_result <- poisson_result[,c(6,5,1:4)]
xtable(poisson_result)
