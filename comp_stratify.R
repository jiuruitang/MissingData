## Created: 11/30
## Modified: 12/01
source("allFunc.R")
doGibbsANWeight<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
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
    Y_mat <- model.matrix(X1~as.factor(Y)+as.factor(S))
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
    #   X1 <- proposed_X1
    #  X1_accepted = X1_accepted + 1
    #}
    
    
    if (i > burnin){
      GAMMA[i-burnin,] <- gamma
      ALPHA[i-burnin,] <- alpha
      X1_impute[i-burnin,] <- X1
    }
  }
  XL <- X1_impute[seq(1,(niter-burnin),100),]
  return(list(gamma = GAMMA,alpha = ALPHA,MI_data = XL, acceptRatio = X1_accepted/niter))
}

alpha_s <- rnorm(4)
gamma_s <- rnorm(4)
testList1 <- doGibbsANWeight()
resultList1 <- getResults(dataMI = testList1$MI_data, n1 = 1500, n2 = 3500,sub_dat)
mean(resultList1$pop_total)
colMeans(resultList1$alpha)
colMeans(resultList1$gamma)

alpha_MAR <- matrix(NA,10,3)
gamma_MAR <- matrix(NA,10,4)
total_MAR <- rep(NA,10)
alpha_ANC <- matrix(NA,10,3)
gamma_ANC <- matrix(NA,10,4)
total_ANC <- rep(NA,10)
acceptRatio_ANC <- rep(NA,10)
alpha_R <- matrix(NA,10,3)
gamma_R <- matrix(NA,10,4)
total_R <- rep(NA,10)
acceptRatio_R <- rep(NA,10)
for(i in 1:10){
  sub_dat <- getSample(population = pop_dat)
  alpha_s <- rnorm(4)
  gamma_s <- rnorm(4)
  testListMAR <- doGibbsMarWeight()
  alpha_s <- rnorm(3)
  testListANC <- doGibbs()
  testListR <- doGibbsSwitch()
  
  resultListMAR <- getResults(dataMI = testListMAR$MI_data, n1 = 1500, n2 = 3500,sub_dat)
  alpha_MAR[i,] <- colMeans(resultListMAR$alpha)
  gamma_MAR[i,] <- colMeans(resultListMAR$gamma)
  total_MAR[i] <- mean(resultListMAR$pop_total)
  
  resultListANC <- getResults(dataMI = testListANC$MI_data, n1 = 1500, n2 = 3500,sub_dat)
  alpha_ANC[i,] <- colMeans(resultListANC$alpha)
  gamma_ANC[i,] <- colMeans(resultListANC$gamma)
  total_ANC[i] <- mean(resultListANC$pop_total)
  acceptRatio_ANC[i] <- testListANC$acceptRatio
  
  resultListR <- getResults(dataMI = testListR$MI_data, n1 = 1500, n2 = 3500,sub_dat)
  alpha_R[i,] <- colMeans(resultListR$alpha)
  gamma_R[i,] <- colMeans(resultListR$gamma)
  total_R[i] <- mean(resultListR$pop_total)
  acceptRatio_R[i] <- testListR$acceptRatio
}
colMeans(alpha_10)
apply(alpha_10,2,sd)/sqrt(10)
colMeans(gamma_10)
apply(gamma_10,2,sd)/sqrt(10)
mean(total_10)
sd(total_10)/sqrt(10)
sum(pop_dat$X1) # truth
mean(acceptRatio_10)

resultMat <- matrix(NA,8,6)
resultMat[,1] <-  c(colMeans(alpha_MAR),colMeans(gamma_MAR),mean(total_MAR))
resultMat[,2] <-  c(apply(alpha_MAR,2,sd)/sqrt(10),apply(gamma_MAR,2,sd)/sqrt(10),sd(total_MAR))
resultMat[,3] <-  c(colMeans(alpha_ANC),colMeans(gamma_ANC),mean(total_ANC))
resultMat[,4] <-  c(apply(alpha_ANC,2,sd)/sqrt(10),apply(gamma_ANC,2,sd)/sqrt(10),sd(total_ANC))
resultMat[,5] <-  c(colMeans(alpha_R),colMeans(gamma_R),mean(total_R))
resultMat[,6] <-  c(apply(alpha_R,2,sd)/sqrt(10),apply(gamma_R,2,sd)/sqrt(10),sd(total_R))


strata_result <- as.data.frame(resultMat)
names(strata_result) <- c("MAR+W Mean","MAR+W SD","AN+C Mean", "AC+C SD","Random+C Mean", "Random+C SD")
strata_result$truth <- c(c(0.5,-0.5,-1),c(-0.25,0.1,0.3,-1.1),sum(pop_dat$X1))
strata_result <- strata_result[,c(7,1:6)]
xtable(strata_result)
