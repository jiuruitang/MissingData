gen_pop_Poisson <- function(N=50000,N1=35000, N2= 15000,alpha = c(0.3,-0.5),gamma = c(-0.25,0.1,-1.1),beta=c(0.5,-0.15)){
  pi <- rbeta(N,2,5) 
  W <- 1/pi
  Y <- W
  
  # generate Y give weights, have some pi_Y very close to 0
  pi_Y <- pnorm(model.matrix(Y~pbeta)%*%beta)
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

pop_dat_Poisson <- gen_pop_Poisson()
getSample_Poisson <- function(population = pop_dat_Poisson,n=5000,N=50000){
  getSampled <- rbinom(N,1,pop_dat_Poisson$pi)
  sub_dat <- pop_dat_Poisson[getSampled==1,]
  return(sub_dat)
}

getSampleStatsPoisson <- function(population = pop_dat_Poisson,N=50000){
  getSampled <- rbinom(N,1,pop_dat_Poisson$pi)
  sub_dat <- pop_dat_Poisson[getSampled==1,]
  return(sum(sub_dat$X1*sub_dat$W))
}
# sample multiple times from the sample
# and get pop mean and variance
sampleTotal <- replicate(200,getSampleStatsPoisson())
pop_sd_HT <- sd(sampleTotal)
pop_mean_HT <- mean(sampleTotal)

sub_dat <- getSample_Poisson()

# sanity check of the parameters in sampled data set
## sanity check passed, similar to truth
m1 <- glm(Rx~as.factor(Y)+as.factor(X1),data=sub_dat,family=binomial(probit))
summary(m1)
sub_dat$id <- row.names(sub_dat)
mydesign <- svydesign(id = ~id,data = sub_dat,weight = ~W)
m2 <- svyglm(X1~as.factor(Y),design = mydesign,family = quasibinomial(link = "probit"))
summary(m2)

alpha_s <- rnorm(2)
beta_s <- rnorm(2)
gamma_s <- rnorm(3)
doGibbs_Poisson<-function(niter = 10000,burnin = 5000, data = sub_dat, alpha = alpha_s,beta = beta_s, gamma = gamma_s, pop_mean = pop_mean_HT, pop_sd = pop_sd_HT){
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
    Z_mean <- alpha[1] + alpha[2]*(Y==1) 
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
    G_mean <- gamma[1] + gamma[2]*Y + gamma[3]*X1
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
    pi_x1 <- pnorm(0,alpha[1] + alpha[2]*Y,1) 
    pi_Rx_temp0 <- dnorm(G,gamma[1] + gamma[2]*Y,1)
    pi_Rx_temp1 <- dnorm(G,gamma[1] + gamma[2]*Y + gamma[3],1)
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

testList <- doGibbs_Poisson()
getResults <- function(dataMI,sub_dat){
  n <- dim(dataMI)[1]
  ans <- matrix(NA,n,2)
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
  ans_g <- matrix(NA,n,3)
  for (i in 1:n){
    test_dat <- as.data.frame(cbind(dataMI[i,],rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$Rx),stringsAsFactors = FALSE)
    names(test_dat) <- c("X1","id","W","Y","Rx")
    test_dat$X1 <- as.factor(test_dat$X1)
    test_dat$id <- as.numeric(test_dat$id)
    test_dat$W <- as.numeric(test_dat$W)
    test_dat$Y <- as.factor(test_dat$Y)
    test_dat$Rx <- as.factor(test_dat$Rx)
    # mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W,strata = ~str)
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

resultList <- getResults(dataMI = testList$MI_data,sub_dat)
colMeans(resultList$alpha)
colMeans(resultList$gamma)
mean(resultList$pop_total)

# repeat this process 10 times (for 10 different population)
alpha_10 <- matrix(NA,10,2)
gamma_10 <- matrix(NA,10,3)
total_10 <- rep(NA,10)
acceptRatio_10 <- rep(NA,10)
pop_dat_Poisson <- gen_pop_Poisson()
for(i in 1:10){
  sampleTotal <- replicate(200,getSampleStatsPoisson())
  pop_sd_HT <- sd(sampleTotal)
  pop_mean_HT <- mean(sampleTotal)
  
  sub_dat <- getSample_Poisson()
  alpha_s <- rnorm(2)
  gamma_s <- rnorm(3)
  testList <- doGibbs_Poisson()
  
  resultList <- getResults(dataMI = testList$MI_data,sub_dat)
  alpha_10[i,] <- colMeans(resultList$alpha)
  gamma_10[i,] <- colMeans(resultList$gamma)
  total_10[i] <- mean(resultList$pop_total)
  acceptRatio_10[i] <- testList$acceptRatio
}
colMeans(alpha_10)
apply(alpha_10,2,sd)/sqrt(10)
colMeans(gamma_10)
apply(gamma_10,2,sd)/sqrt(10)
mean(total_10)
sd(total_10)/sqrt(10)
sum(pop_dat_Poisson$X1) # truth
mean(acceptRatio_10)

boxplot(log(W)~Y,data = pop_dat_Poisson)
pop_dat_Poisson %>% group_by(X1) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
