gen_pop_Poisson <- function(N=50000,N1=35000, N2= 15000,alpha = c(0.3,-0.5),gamma = c(-0.25,0.1,-1.1),beta=c(0.5,-0.15)){
  pbeta <- rbeta(N,2,5) 
  W <- 1/pbeta
  Y <- W
  
  # beta=c(0.5,2.8)
  # alpha = c(0.3,-0.5,2.5)
  gamma = c(-1.3,0.1,-2.7, 3.6)
  # generate Y give weights, have some pi_Y very close to 0
  pi_Y <- pnorm(model.matrix(Y~pbeta)%*%beta)
  Y <- rbinom(N,1,p=pi_Y)
  
  # generate X1 given Y
  X1 <- Y
  pi_x1 <- pnorm(model.matrix(X1~Y)%*%alpha) # pnorm is CDF or normal
  X1 <- rbinom(N,1,p=pi_x1)
  
  # sample Rx|X,Y
  Rx <- X1
  pi_Rx <- pnorm(model.matrix(Rx~Y+X1+pbeta)%*%gamma)
  Rx <- rbinom(N,1,p = pi_Rx)
  
  # get a dataframe for population
  return(as.data.frame(cbind(Y,X1,Rx,W,pbeta)))
}

pop_dat_Poisson <- gen_pop_Poisson()
mean(pop_dat_Poisson$Rx)
# boxplot(log(W)~Y,data = pop_dat_Poisson)
library(xtable)
YW <- pop_dat_Poisson %>% group_by(Y) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
X1W <- pop_dat_Poisson %>% group_by(X1) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
RxW <- pop_dat_Poisson %>% group_by(Rx) %>% summarise_at(vars(W), list(~mean(., na.rm=TRUE)))
YW;X1W;RxW
xtable(YW)
xtable(X1W)
xtable(RxW)
table(pop_dat_Poisson$Y,pop_dat_Poisson$X1)
summary(glm(X1~Y,data = pop_dat_Poisson,family = "binomial"))
