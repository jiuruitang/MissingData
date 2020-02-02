## Created: 12/01
## Try to compare three methods with different Poisson sampling setting

source("allFunc.R")

## Setting 1: Y depends on 1/weight, Rx and X1 depends of Y
pop_dat <- gen_pop_Poisson()
sub_dat <- getSample_Poisson(population = pop_dat)
sampleTotal <- replicate(200,getSampleStatsPoisson(population = pop_dat))
pop_sd_HT <- sd(sampleTotal)
pop_mean_HT <- mean(sampleTotal)

alpha_s <- rnorm(3)
gamma_s <- rnorm(3)
testListMAR <- doGibbsMarWeightPoisson()
resultListMAR <- getResults(dataMI = testListMAR$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
colMeans(resultListMAR$alpha)
colMeans(resultListMAR$gamma)
mean(resultListMAR$pop_total)
sum(pop_dat$X1) # truth

testListANWC <- doGibbsANWCPoisson()
resultListANWC <- getResults(dataMI = testListANWC$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
colMeans(resultListANWC$alpha)
colMeans(resultListANWC$gamma)
mean(resultListANWC$pop_total)

testListSwitch <- doGibbsSwitchWCPoisson()
resultListSwitch <- getResults(dataMI = testListSwitch$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
colMeans(resultListSwitch$alpha)
colMeans(resultListSwitch$gamma)
mean(resultListSwitch$pop_total)

mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
m1 <- svyglm(X1~Y+ I(1/W),design = mydesign,family = quasibinomial(link = "probit"))
coef(m1)
m2 <- glm(Rx~Y+X1+I(1/W),data = sub_dat, family = binomial(link = "probit"))


alpha_MAR <- matrix(NA,10,2)
gamma_MAR <- matrix(NA,10,3)
total_MAR <- rep(NA,10)
alpha_ANWC <- matrix(NA,10,2)
gamma_ANWC <- matrix(NA,10,3)
total_ANWC <- rep(NA,10)
acceptRatio_ANWC <- rep(NA,10)
alpha_R <- matrix(NA,10,2)
gamma_R <- matrix(NA,10,3)
total_R <- rep(NA,10)
acceptRatio_R <- rep(NA,10)
alpha_truth <- matrix(NA,10,2)
gamma_truth <- matrix(NA,10,3)
for(i in 2:10){
  sub_dat <- getSample_Poisson(population = pop_dat)
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
  
  testListR <- doGibbsSwitchWCPoisson()
  resultListR <- getResults(dataMI = testListR$MI_data,alpha_len=2, gamma_len = 3,sub_dat)
  alpha_R[i,] <- colMeans(resultListR$alpha)
  gamma_R[i,] <- colMeans(resultListR$gamma)
  total_R[i] <- mean(resultListR$pop_total)
  acceptRatio_R[i] <- testListR$acceptRatio
}

resultMat <- matrix(NA,6,6)
resultMat[,1] <-  c(colMeans(alpha_MAR),colMeans(gamma_MAR),mean(total_MAR))
resultMat[,2] <-  c(apply(alpha_MAR,2,sd)/sqrt(10),apply(gamma_MAR,2,sd)/sqrt(10),sd(total_MAR))
resultMat[,3] <-  c(colMeans(alpha_ANWC),colMeans(gamma_ANWC),mean(total_ANWC))
resultMat[,4] <-  c(apply(alpha_ANWC,2,sd)/sqrt(10),apply(gamma_ANWC,2,sd)/sqrt(10),sd(total_ANWC))
resultMat[,5] <-  c(colMeans(alpha_R),colMeans(gamma_R),mean(total_R))
resultMat[,6] <-  c(apply(alpha_R,2,sd)/sqrt(10),apply(gamma_R,2,sd)/sqrt(10),sd(total_R))

poisson_result <- as.data.frame(resultMat)
names(poisson_result) <- c("MAR+W Mean","MAR+W SD","AN+WC Mean", "AN+WC SD","Random+C Mean", "Random+C SD")
poisson_result$Sampletruth <- c(colMeans(alpha_truth),colMeans(gamma_truth),sum(pop_dat$X1))
poisson_result$truth <-c(c(0.3,-0.5),c(-0.25,0.1,-1.1),sum(pop_dat$X1))
poisson_result <- poisson_result[,c(8,7,1:6)]
xtable(poisson_result)

## Setting 2: Weight depends on X and Y
pop_dat <- gen_pop_PPS()

sampleTotal <- replicate(200,getSampleStatsPoisson(population = pop_dat))
pop_sd_HT <- sd(sampleTotal)
pop_mean_HT <- mean(sampleTotal)

alpha_MARP <- matrix(NA,10,3)
gamma_MARP <- matrix(NA,10,4)
total_MARP <- rep(NA,10)
alpha_ANWCP <- matrix(NA,10,3)
gamma_ANWCP <- matrix(NA,10,4)
total_ANWCP <- rep(NA,10)
acceptRatio_ANWCP <- rep(NA,10)
alpha_RP <- matrix(NA,10,3)
gamma_RP <- matrix(NA,10,4)
total_RP <- rep(NA,10)
total_truth <- rep(NA,10)
acceptRatio_RP <- rep(NA,10)
alpha_truthP <- matrix(NA,10,3)
gamma_truthP <- matrix(NA,10,4)
for(i in 1:10){
  sub_dat <- getSample_Poisson(population = pop_dat)
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y+ I(1/W),design = mydesign,family = quasibinomial(link = "probit"))
  alpha_truthP[i,] <- coef(m1)[1:3]
  m2 <- glm(Rx~Y+X1+I(1/W),data = sub_dat, family = binomial(link = "probit"))
  gamma_truthP[i,] <- coef(m2)[1:4]
  total_truth[i] <- sum(sub_dat$X1*sub_dat$W)
  
  alpha_s <- rnorm(4)
  gamma_s <- rnorm(4)
  testListMAR <- doGibbsMarWeightPoisson()
  resultListMAR <- getResults(dataMI = testListMAR$MI_data,alpha_len=3, gamma_len = 4,sub_dat)
  alpha_MARP[i,] <- colMeans(resultListMAR$alpha)
  gamma_MARP[i,] <- colMeans(resultListMAR$gamma)
  total_MARP[i] <- mean(resultListMAR$pop_total)
  
  testListANWC <- doGibbsANWCPoisson()
  resultListANWC <- getResults(dataMI = testListANWC$MI_data,alpha_len=3, gamma_len = 4,sub_dat)
  alpha_ANWCP[i,] <- colMeans(resultListANWC$alpha)
  gamma_ANWCP[i,] <- colMeans(resultListANWC$gamma)
  total_ANWCP[i] <- mean(resultListANWC$pop_total)
  acceptRatio_ANWCP[i] <- testListANWC$acceptRatio
  
  testListR <- doGibbsSwitchWCPoisson()
  resultListR <- getResults(dataMI = testListR$MI_data,alpha_len=3, gamma_len = 4, sub_dat)
  alpha_RP[i,] <- colMeans(resultListR$alpha)
  gamma_RP[i,] <- colMeans(resultListR$gamma)
  total_RP[i] <- mean(resultListR$pop_total)
  acceptRatio_RP[i] <- testListR$acceptRatio
}

resultMatP <- matrix(NA,8,6)
resultMatP[,1] <-  c(colMeans(alpha_MARP),colMeans(gamma_MARP),mean(total_MARP))
resultMatP[,2] <-  c(apply(alpha_MARP,2,sd)/sqrt(10),apply(gamma_MARP,2,sd)/sqrt(10),sd(total_MARP))
resultMatP[,3] <-  c(colMeans(alpha_ANWCP),colMeans(gamma_ANWCP),mean(total_ANWCP))
resultMatP[,4] <-  c(apply(alpha_ANWCP,2,sd)/sqrt(10),apply(gamma_ANWCP,2,sd)/sqrt(10),sd(total_ANWCP))
resultMatP[,5] <-  c(colMeans(alpha_RP),colMeans(gamma_RP),mean(total_RP))
resultMatP[,6] <-  c(apply(alpha_RP,2,sd)/sqrt(10),apply(gamma_RP,2,sd)/sqrt(10),sd(total_RP))

PPS_result <- as.data.frame(resultMatP)
names(PPS_result) <- c("MAR+W Mean","MAR+W SD","AN+WC Mean", "AN+WC SD","Random+C Mean", "Random+C SD")
PPS_result$Sampletruth <- c(colMeans(alpha_truthP),colMeans(gamma_truthP),mean(total_truth))
PPS_result$truth <-c(c(0.5,-0.5,-1),c(-0.6,0.1,0.3,-1.1),sum(pop_dat$X1))
PPS_result <- PPS_result[,c(8,7,1:6)]
PPS_resultR <- round(PPS_result,digits=5)
xtable(PPS_result,digits = -2)
