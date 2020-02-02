source("allFunc.R")
gen_pop <- function(N=50000,alpha = c(0.5,-0.5,-1),gamma = c(-0.6,0.1,0.3,-1.1),beta= c(2,-1,0.8),theta = c(0.1,0.45,0.45)){
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


pop_dat <- gen_pop()

YW <- pop_dat %>% group_by(Y) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
X1W <- pop_dat %>% group_by(X1) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
RxW <- pop_dat %>% group_by(Rx) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
YW;X1W;RxW
xtable(YW)
xtable(X1W)
xtable(RxW)
sampleTotal <- replicate(200,getSampleStatsPoisson(population = pop_dat))
pop_sd_HT <- sd(sampleTotal)
pop_mean_HT <- mean(sampleTotal)

alpha_10 <- matrix(NA,10,3)
gamma_10 <- matrix(NA,10,4)
total_10 <- rep(NA,10)
acceptRatio_10 <- rep(NA,10)
for(i in 1:10){
  sub_dat <- getSample(population = pop_dat)
  alpha_s <- rnorm(4)
  gamma_s <- rnorm(4)
  testList <- doGibbsANWeight()
  
  resultList <- getResults(dataMI = testList$MI_data, n1 = 1500, n2 = 3500,sub_dat)
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
sum(pop_dat$X1) # truth
mean(acceptRatio_10)

# try random switch
alpha_10_2 <- matrix(NA,10,3)
gamma_10_2 <- matrix(NA,10,4)
total_10_2 <- rep(NA,10)
acceptRatio_10_2 <- rep(NA,10)
for(i in 1:10){
  sub_dat <- getSample_Poisson(population = pop_dat)
  alpha_s <- rnorm(3)
  gamma_s <- rnorm(4)
  testList <- doGibbsSwitchGoodStart()
  
  resultList <- getResults(dataMI = testList$MI_data, n1 = 1500, n2 = 3500,sub_dat)
  alpha_10_2[i,] <- colMeans(resultList$alpha)
  gamma_10_2[i,] <- colMeans(resultList$gamma)
  total_10_2[i] <- mean(resultList$pop_total)
  acceptRatio_10_2[i] <- testList$acceptRatio
}
colMeans(alpha_10_2)
apply(alpha_10_2,2,sd)/sqrt(10)
colMeans(gamma_10_2)
apply(gamma_10_2,2,sd)/sqrt(10)
mean(total_10_2)
sd(total_10_2)/sqrt(10)
sum(pop_dat$X1) # truth
mean(acceptRatio_10_2)

testList_M <- doGibbs()
plot(testList$gamma[4900:5000,1],type="l")
plot(testList_M$gamma[4900:5000,1],type="l")
plot(testList$gamma[,2],type="l")
plot(testList_M$gamma[,2],type="l")

sub_dat$W[which(testList$MI_data[1,][which(Rx == 1)] != sub_dat$X1[which(Rx == 1)])]
sub_dat$X1[which(testList$MI_data[1,][which(Rx == 1)] != sub_dat$X1[which(Rx == 1)])]

junk_var <- sum((pop_dat$Y/pop_dat$pi)^2*pop_dat$pi*(1-pop_dat$pi))
