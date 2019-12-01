source("allFunc.R")
pop_dat <- gen_pop()
sampleTotal <- replicate(200,getSampleStats(population = pop_dat))
pop_sd_HT <- sd(sampleTotal)
pop_mean_HT <- mean(sampleTotal)

sub_dat <- getSample(population = pop_dat)

pop_dat$pi <- 1/pop_dat$W
sub_dat$pi <- 1/sub_dat$W

# calculate variance based on HT-estimator
pi_ij <- matrix(as.numeric(outer(sub_dat$W,sub_dat$W, FUN= "==")),nrow=5000)
piipij <- matrix(outer(sub_dat$pi,sub_dat$pi,FUN="*"),nrow=5000)
yiyj <- matrix(outer(sub_dat$X1,sub_dat$X1,FUN="*"),nrow=5000)
junk_var <- sum((sub_dat$X1/sub_dat$pi)^2*sub_dat$pi*(1-sub_dat$pi)) + sum(yiyj/piipij*(pi_ij - piipij))

# start at truth for stratified sampling
alpha_10 <- matrix(NA,10,3)
gamma_10 <- matrix(NA,10,4)
total_10 <- rep(NA,10)
acceptRatio_10 <- rep(NA,10)
alpha_10_2 <- matrix(NA,10,3)
gamma_10_2 <- matrix(NA,10,4)
total_10_2 <- rep(NA,10)
acceptRatio_10_2 <- rep(NA,10)
for(i in 1:10){
  sub_dat <- getSample(population = pop_dat)
  alpha_s <- rnorm(3)
  gamma_s <- rnorm(4)
  testList1 <- doGibbsGoodStart()
  
  resultList1 <- getResults(dataMI = testList1$MI_data, n1 = 1500, n2 = 3500,sub_dat)
  alpha_10[i,] <- colMeans(resultList1$alpha)
  gamma_10[i,] <- colMeans(resultList1$gamma)
  total_10[i] <- mean(resultList1$pop_total)
  acceptRatio_10[i] <- testList1$acceptRatio
  
  testList2 <- doGibbsSwitchGoodStart()
  
  resultList2 <- getResults(dataMI = testList2$MI_data, n1 = 1500, n2 = 3500,sub_dat)
  alpha_10_2[i,] <- colMeans(resultList2$alpha)
  gamma_10_2[i,] <- colMeans(resultList2$gamma)
  total_10_2[i] <- mean(resultList2$pop_total)
  acceptRatio_10_2[i] <- testList2$acceptRatio
  
}
colMeans(alpha_10)
apply(alpha_10,2,sd)/sqrt(10)
colMeans(gamma_10)
apply(gamma_10,2,sd)/sqrt(10)
mean(total_10)
sd(total_10)/sqrt(10)
sum(pop_dat$X1) # truth
mean(acceptRatio_10)

colMeans(alpha_10_2)
apply(alpha_10_2,2,sd)/sqrt(10)
colMeans(gamma_10_2)
apply(gamma_10_2,2,sd)/sqrt(10)
mean(total_10_2)
sd(total_10_2)/sqrt(10)
sum(pop_dat$X1) # truth
mean(acceptRatio_10_2)
