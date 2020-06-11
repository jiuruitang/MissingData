# This script analyzes the results from CPS_server.R
# created: 02/03
# modified: 02/19
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss2 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
for (i in 1:sim_n){
  load(paste("~/Missing/results2/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- mean(resultListANWC$total_mean)
  total <- resultListANWC$total_mean
  total_var <- resultListANWC$total_var
  u_L_t <- mean(total_var)
  b_L_t <- apply((scale(total,scale=FALSE))^2/(n-1),2,sum)
  ANWC_total_var[i] <- (1+1/n)*b_L_t + u_L_t
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss2[i,] <- coef(m1)
}

library(plotrix)
plotCI(x=1:sim_n,y=ANWC_mean[,1],li=ANWC_mean[,1]-2.009*sqrt(ANWC_var[,1]),ui=ANWC_mean[,1]+2.009*sqrt(ANWC_var[,1]),main = expression(alpha[0])) # plot the confidence intervals. The first row of ci, namely ci[1,], has the left end points. The second row of ci, namely ci[2,], has the right end points.
abline(h=0.3,col="red") # draw a horizontal line at the value of the true mean
sum(ANWC_mean[,1]-2.009*sqrt(ANWC_var[,1]) > 0.3 | ANWC_mean[,1]+2.009*sqrt(ANWC_var[,1]) < 0.3)

plotCI(x=1:sim_n,y=ANWC_mean[,2],li=ANWC_mean[,2]-2.009*sqrt(ANWC_var[,2]),ui=ANWC_mean[,2]+2.009*sqrt(ANWC_var[,2]),main = expression(alpha[1])) # plot the confidence intervals. The first row of ci, namely ci[1,], has the left end points. The second row of ci, namely ci[2,], has the right end points.
abline(h=-0.5,col="red") # draw a horizontal line at the value of the true mean
sum(ANWC_mean[,2]-2.009*sqrt(ANWC_var[,2]) > -0.5 | ANWC_mean[,2]+2.009*sqrt(ANWC_var[,2]) < -0.5)

plotCI(x=1:sim_n,y=ANWC_mean[,3],li=ANWC_mean[,3]-2.009*sqrt(ANWC_var[,3]),ui=ANWC_mean[,3]+2.009*sqrt(ANWC_var[,3]),main = expression(gamma[0])) # plot the confidence intervals. The first row of ci, namely ci[1,], has the left end points. The second row of ci, namely ci[2,], has the right end points.
abline(h=-0.25,col="red") # draw a horizontal line at the value of the true mean
sum(ANWC_mean[,3]-2.009*sqrt(ANWC_var[,3]) > -0.25 | ANWC_mean[,3]+2.009*sqrt(ANWC_var[,3]) < -0.25)

plotCI(x=1:sim_n,y=ANWC_mean[,4],li=ANWC_mean[,4]-2.009*sqrt(ANWC_var[,4]),ui=ANWC_mean[,4]+2.009*sqrt(ANWC_var[,4]),main = expression(gamma[1])) # plot the confidence intervals. The first row of ci, namely ci[1,], has the left end points. The second row of ci, namely ci[2,], has the right end points.
abline(h=0.1,col="red") # draw a horizontal line at the value of the true mean
sum(ANWC_mean[,4]-2.009*sqrt(ANWC_var[,4]) > 0.1 | ANWC_mean[,4]+2.009*sqrt(ANWC_var[,4]) < 0.1)

plotCI(x=1:sim_n,y=ANWC_mean[,5],li=ANWC_mean[,5]-2.009*sqrt(ANWC_var[,5]),ui=ANWC_mean[,5]+2.009*sqrt(ANWC_var[,5]),main = expression(gamma[2])) # plot the confidence intervals. The first row of ci, namely ci[1,], has the left end points. The second row of ci, namely ci[2,], has the right end points.
abline(h=-1.5,col="red") # draw a horizontal line at the value of the true mean
sum(ANWC_mean[,5]-2.009*sqrt(ANWC_var[,5]) > -1.5 | ANWC_mean[,5]+2.009*sqrt(ANWC_var[,5]) < -1.5)


# coverage of CI
res_df2 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df2$ANWC_mean.1 - 2.009*sqrt(res_df2$ANWC_var.1) > 0.3|res_df2$ANWC_mean.1 + 2.009*sqrt(res_df2$ANWC_var.1) < 0.3)/1000
1-sum(res_df2$ANWC_mean.2 - 2.009*sqrt(res_df2$ANWC_var.2) > -0.5|res_df2$ANWC_mean.2 + 2.009*sqrt(res_df2$ANWC_var.2) < -0.5)/1000
1-sum(res_df2$ANWC_mean.3 - 2.009*sqrt(res_df2$ANWC_var.3) > -0.25|res_df2$ANWC_mean.3 + 2.009*sqrt(res_df2$ANWC_var.3) < -0.25)/1000
1-sum(res_df2$ANWC_mean.4 - 2.009*sqrt(res_df2$ANWC_var.4) > 0.1|res_df2$ANWC_mean.4 + 2.009*sqrt(res_df2$ANWC_var.4) < 0.1)/1000
1-sum(res_df2$ANWC_mean.5 - 2.009*sqrt(res_df2$ANWC_var.5) > -1.5|res_df2$ANWC_mean.5 + 2.009*sqrt(res_df2$ANWC_var.5) < -1.5)/1000
1-sum(res_df2$ANWC_T - 2.009*sqrt(res_df2$ANWC_TV) > 64110 | res_df2$ANWC_T + 2.009*sqrt(res_df2$ANWC_TV) < 64110)/1000

# MCMC variance
mean(res_df$ANWC_var.5)
var(res_df$ANWC_mean.5)
mean(res_df$ANWC_var.4)
var(res_df$ANWC_mean.4)
mean(res_df$ANWC_var.3)
var(res_df$ANWC_mean.3)
mean(res_df$ANWC_var.2)
var(res_df$ANWC_mean.2)
mean(res_df$ANWC_var.1)
var(res_df$ANWC_mean.1)
mean(res_df$ANWC_TV)
var(res_df$ANWC_T)

########################
# 02/16 Update
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss3 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
for (i in 1:sim_n){
  load(paste("~/Missing/results3/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss3[i,] <- coef(m1)
}

# coverage of CI
res_df3 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df3$ANWC_mean.1 - 2.009*sqrt(res_df3$ANWC_var.1) > 0.3|res_df3$ANWC_mean.1 + 2.009*sqrt(res_df3$ANWC_var.1) < 0.3)/1000
1-sum(res_df3$ANWC_mean.2 - 2.009*sqrt(res_df3$ANWC_var.2) > -0.5|res_df3$ANWC_mean.2 + 2.009*sqrt(res_df3$ANWC_var.2) < -0.5)/1000
1-sum(res_df3$ANWC_mean.3 - 2.009*sqrt(res_df3$ANWC_var.3) > -0.25|res_df3$ANWC_mean.3 + 2.009*sqrt(res_df3$ANWC_var.3) < -0.25)/1000
1-sum(res_df3$ANWC_mean.4 - 2.009*sqrt(res_df3$ANWC_var.4) > 0.1|res_df3$ANWC_mean.4 + 2.009*sqrt(res_df3$ANWC_var.4) < 0.1)/1000
1-sum(res_df3$ANWC_mean.5 - 2.009*sqrt(res_df3$ANWC_var.5) > -1.5|res_df3$ANWC_mean.5 + 2.009*sqrt(res_df3$ANWC_var.5) < -1.5)/1000
1-sum(res_df3$ANWC_T - 2.009*sqrt(res_df3$ANWC_TV) > 64110 | res_df3$ANWC_T + 2.009*sqrt(res_df3$ANWC_TV) < 64110)/1000

# MCMC variance
mean(res_df$ANWC_var.5)
var(res_df$ANWC_mean.5)
mean(res_df$ANWC_var.4)
var(res_df$ANWC_mean.4)
mean(res_df$ANWC_var.3)
var(res_df$ANWC_mean.3)
mean(res_df$ANWC_var.2)
var(res_df$ANWC_mean.2)
mean(res_df$ANWC_var.1)
var(res_df$ANWC_mean.1)
mean(res_df$ANWC_TV)
var(res_df$ANWC_T)

### truly update
which(res_df2$ANWC_var.5 == max(res_df2$ANWC_var.5)) # 411
m2 <- glm(Rx~Y+X1,data = sub_dat, family = binomial(link = "probit"))
mean(resultListANWC$total_mean) # 62161.24

obs_dat <- sub_dat[which(sub_dat$Rx == 0),]
mis_dat <- sub_dat[which(sub_dat$Rx == 1),]
sum(obs_dat$W*obs_dat$X1)
sum(mis_dat$W*mis_dat$X1)
sum(sub_dat$X1*sub_dat$W) # truth : 64929.0

sort(res_df2$ANWC_var.5, decreasing = T)[1:15]
sample_total <- rep(NA,15)
for (i in 1:15){
  var_val <- sort(res_df2$ANWC_var.5, decreasing = T)[i]
  index <- which(res_df2$ANWC_var.5 == var_val)
  load(paste("~/Missing/results2/Mis_",index,".RData",sep=""))
  sample_total[i] <- sum(sub_dat$X1*sub_dat$W)
}
sample_total > (64110 + 3*35) # 13/15 is true

# smaller variance 
which(res_df3$ANWC_var.5 == max(res_df3$ANWC_var.5)) # 575
m2 <- glm(Rx~Y+X1,data = sub_dat, family = binomial(link = "probit"))
summary(m2)

mean(resultListANWC$total_mean) # estimate : 64428.23
obs_dat <- sub_dat[which(sub_dat$Rx == 0),]
mis_dat <- sub_dat[which(sub_dat$Rx == 1),]
sum(obs_dat$W*obs_dat$X1)
sum(mis_dat$W*mis_dat$X1)
sum(sub_dat$X1*sub_dat$W) # sample truth : 67183.59

sample_total <- rep(NA,15)
for (i in 1:15){
  var_val <- sort(res_df3$ANWC_var.5, decreasing = T)[i]
  index <- which(res_df3$ANWC_var.5 == var_val)
  load(paste("~/Missing/results3/Mis_",index,".RData",sep=""))
  sample_total[i] <- sum(sub_dat$X1*sub_dat$W)
}
sample_total > (64110 + 3*35) # all of them is true

test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
names(test_dat) <- c("X1","id","W","Y")
test_dat$X1 <- as.factor(test_dat$X1)
test_dat$id <- as.numeric(test_dat$id)
test_dat$W <- as.numeric(test_dat$W)
test_dat$Y <- as.factor(test_dat$Y)
mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))

# 02/20 Update
load("~/Missing/results3/Mis_575.RData")
n <- dim(MI_dataANWC)[1]
ans_g <- matrix(NA,n,3)
ul_g <- matrix(NA,n,3)

for (i in 1:n){
  test_dat <- as.data.frame(cbind(MI_dataANWC[i,],rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- glm(Rx~Y+X1,data = test_dat, family = binomial(link = "probit"))
  ans_g[i,] <- coef(m2)
  ul_g[i,] <- diag(vcov(m2))
}
# colMeans(ans_g)
u_L_g <- colMeans(ul_g)
b_L_g <- apply((scale(ans_g,scale=FALSE))^2/(n-1),2,sum)
T_L_g <- (1+1/n)*b_L_g + u_L_g

apply(MI_dataANWC[,which(sub_dat$Rx == 1)],1,sum)


# 02/23 Update
# Result of even smaller Vx
# coverage of CI
res_df4 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df4$ANWC_mean.1 - 2.009*sqrt(res_df4$ANWC_var.1) > 0.3|res_df4$ANWC_mean.1 + 2.009*sqrt(res_df4$ANWC_var.1) < 0.3)/1000
1-sum(res_df4$ANWC_mean.2 - 2.009*sqrt(res_df4$ANWC_var.2) > -0.5|res_df4$ANWC_mean.2 + 2.009*sqrt(res_df4$ANWC_var.2) < -0.5)/1000
1-sum(res_df4$ANWC_mean.3 - 2.009*sqrt(res_df4$ANWC_var.3) > -0.25|res_df4$ANWC_mean.3 + 2.009*sqrt(res_df4$ANWC_var.3) < -0.25)/1000
1-sum(res_df4$ANWC_mean.4 - 2.009*sqrt(res_df4$ANWC_var.4) > 0.1|res_df4$ANWC_mean.4 + 2.009*sqrt(res_df4$ANWC_var.4) < 0.1)/1000
1-sum(res_df4$ANWC_mean.5 - 2.009*sqrt(res_df4$ANWC_var.5) > -1.5|res_df4$ANWC_mean.5 + 2.009*sqrt(res_df4$ANWC_var.5) < -1.5)/1000
1-sum(res_df4$ANWC_T - 2.009*sqrt(res_df4$ANWC_TV) > 64110 | res_df4$ANWC_T + 2.009*sqrt(res_df4$ANWC_TV) < 64110)/1000

# MCMC variance
mean(res_df4$ANWC_var.5)
var(res_df4$ANWC_mean.5)
mean(res_df4$ANWC_var.4)
var(res_df4$ANWC_mean.4)
mean(res_df4$ANWC_var.3)
var(res_df4$ANWC_mean.3)
mean(res_df4$ANWC_var.2)
var(res_df4$ANWC_mean.2)
mean(res_df4$ANWC_var.1)
var(res_df4$ANWC_mean.1)
mean(res_df4$ANWC_TV)
var(res_df4$ANWC_T)

sample_total <- rep(NA,15)
for (i in 1:15){
  var_val <- sort(res_df4$ANWC_var.5, decreasing = T)[i]
  index <- which(res_df4$ANWC_var.5 == var_val)
  load(paste("~/Missing/results4/Mis_",index,".RData",sep=""))
  sample_total[i] <- sum(sub_dat$X1*sub_dat$W)
}
sample_total > (64110 + 3*35) # all of them is true

# draw the 45 degree scatter plot of pre-missing and imputation estimate of parameters
alpha_premiss <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
for (i in 1:sim_n){
  load(paste("~/Missing/results4/Mis_",i,".RData",sep=""))
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss[i,] <- coef(m1)
}

plot(res_df4$ANWC_mean.1,alpha_premiss[,1],xlim = c(0.1,0.5),ylim = c(0.1,0.5),pch=16)
abline(0,1,col = "red")
p <- qplot(res_df4$ANWC_mean.1,alpha_premiss[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df4$ANWC_mean.2,alpha_premiss[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df2$ANWC_mean.1,alpha_premiss2[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df2$ANWC_mean.2,alpha_premiss2[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df3$ANWC_mean.1,alpha_premiss3[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df3$ANWC_mean.2,alpha_premiss3[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

############################
###### 03/01 Update   
# process results5
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss5 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results5/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss5[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# coverage of CI
res_df5 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df5$ANWC_mean.1 - 2.009*sqrt(res_df5$ANWC_var.1) > 0.3|res_df5$ANWC_mean.1 + 2.009*sqrt(res_df5$ANWC_var.1) < 0.3)/1000
1-sum(res_df5$ANWC_mean.2 - 2.009*sqrt(res_df5$ANWC_var.2) > -1|res_df5$ANWC_mean.2 + 2.009*sqrt(res_df5$ANWC_var.2) < -1)/1000
1-sum(res_df5$ANWC_mean.3 - 2.009*sqrt(res_df5$ANWC_var.3) > -0.25|res_df5$ANWC_mean.3 + 2.009*sqrt(res_df5$ANWC_var.3) < -0.25)/1000
1-sum(res_df5$ANWC_mean.4 - 2.009*sqrt(res_df5$ANWC_var.4) > 0.1|res_df5$ANWC_mean.4 + 2.009*sqrt(res_df5$ANWC_var.4) < 0.1)/1000
1-sum(res_df5$ANWC_mean.5 - 2.009*sqrt(res_df5$ANWC_var.5) > -1.5|res_df5$ANWC_mean.5 + 2.009*sqrt(res_df5$ANWC_var.5) < -1.5)/1000
1-sum(res_df5$ANWC_T - 2.009*sqrt(res_df5$ANWC_TV) > 48030 | res_df5$ANWC_T + 2.009*sqrt(res_df5$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df5$ANWC_var.5)
var(res_df5$ANWC_mean.5)
mean(res_df5$ANWC_var.4)
var(res_df5$ANWC_mean.4)
mean(res_df5$ANWC_var.3)
var(res_df5$ANWC_mean.3)
mean(res_df5$ANWC_var.2)
var(res_df5$ANWC_mean.2)
mean(res_df5$ANWC_var.1)
var(res_df5$ANWC_mean.1)
mean(res_df5$ANWC_TV)
var(res_df5$ANWC_T)

p <- qplot(res_df5$ANWC_mean.1,alpha_premiss5[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df5$ANWC_mean.2,alpha_premiss5[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

# process results6
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss6 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results6/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss6[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# coverage of CI
res_df6 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df6$ANWC_mean.1 - 2.009*sqrt(res_df6$ANWC_var.1) > 0.3|res_df6$ANWC_mean.1 + 2.009*sqrt(res_df6$ANWC_var.1) < 0.3)/1000
1-sum(res_df6$ANWC_mean.2 - 2.009*sqrt(res_df6$ANWC_var.2) > -1|res_df6$ANWC_mean.2 + 2.009*sqrt(res_df6$ANWC_var.2) < -1)/1000
1-sum(res_df6$ANWC_mean.3 - 2.009*sqrt(res_df6$ANWC_var.3) > -0.25|res_df6$ANWC_mean.3 + 2.009*sqrt(res_df6$ANWC_var.3) < -0.25)/1000
1-sum(res_df6$ANWC_mean.4 - 2.009*sqrt(res_df6$ANWC_var.4) > 0.1|res_df6$ANWC_mean.4 + 2.009*sqrt(res_df6$ANWC_var.4) < 0.1)/1000
1-sum(res_df6$ANWC_mean.5 - 2.009*sqrt(res_df6$ANWC_var.5) > -1.5|res_df6$ANWC_mean.5 + 2.009*sqrt(res_df6$ANWC_var.5) < -1.5)/1000
1-sum(res_df6$ANWC_T - 2.009*sqrt(res_df6$ANWC_TV) > 48030 | res_df6$ANWC_T + 2.009*sqrt(res_df6$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df6$ANWC_var.5)
var(res_df6$ANWC_mean.5)
mean(res_df6$ANWC_var.4)
var(res_df6$ANWC_mean.4)
mean(res_df6$ANWC_var.3)
var(res_df6$ANWC_mean.3)
mean(res_df6$ANWC_var.2)
var(res_df6$ANWC_mean.2)
mean(res_df6$ANWC_var.1)
var(res_df6$ANWC_mean.1)
mean(res_df6$ANWC_TV)
var(res_df6$ANWC_T)

p <- qplot(res_df6$ANWC_mean.1,alpha_premiss6[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df6$ANWC_mean.2,alpha_premiss6[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100


# process results7
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss7 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results7/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss7[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# coverage of CI
res_df7 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df7$ANWC_mean.1 - 2.009*sqrt(res_df7$ANWC_var.1) > 0.3|res_df7$ANWC_mean.1 + 2.009*sqrt(res_df7$ANWC_var.1) < 0.3)/1000
1-sum(res_df7$ANWC_mean.2 - 2.009*sqrt(res_df7$ANWC_var.2) > -1|res_df7$ANWC_mean.2 + 2.009*sqrt(res_df7$ANWC_var.2) < -1)/1000
1-sum(res_df7$ANWC_mean.3 - 2.009*sqrt(res_df7$ANWC_var.3) > -0.25|res_df7$ANWC_mean.3 + 2.009*sqrt(res_df7$ANWC_var.3) < -0.25)/1000
1-sum(res_df7$ANWC_mean.4 - 2.009*sqrt(res_df7$ANWC_var.4) > 0.1|res_df7$ANWC_mean.4 + 2.009*sqrt(res_df7$ANWC_var.4) < 0.1)/1000
1-sum(res_df7$ANWC_mean.5 - 2.009*sqrt(res_df7$ANWC_var.5) > -1.5|res_df7$ANWC_mean.5 + 2.009*sqrt(res_df7$ANWC_var.5) < -1.5)/1000
1-sum(res_df7$ANWC_T - 2.009*sqrt(res_df7$ANWC_TV) > 48030 | res_df7$ANWC_T + 2.009*sqrt(res_df7$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df7$ANWC_var.5)
var(res_df7$ANWC_mean.5)
mean(res_df7$ANWC_var.4)
var(res_df7$ANWC_mean.4)
mean(res_df7$ANWC_var.3)
var(res_df7$ANWC_mean.3)
mean(res_df7$ANWC_var.2)
var(res_df7$ANWC_mean.2)
mean(res_df7$ANWC_var.1)
var(res_df7$ANWC_mean.1)
mean(res_df7$ANWC_TV)
var(res_df7$ANWC_T)

p <- qplot(res_df7$ANWC_mean.1,alpha_premiss7[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df7$ANWC_mean.2,alpha_premiss7[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100

save(res_df7,alpha_premiss7,TPR,TNR,FPR,FNR,file="./results7/res_df7.RData")


TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results4/Mis_",i,".RData",sep=""))
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# 03/03 Update
# process results8
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss8 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results8/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss8[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# coverage of CI
res_df8 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df8$ANWC_mean.1 - 2.009*sqrt(res_df8$ANWC_var.1) > 0.3|res_df8$ANWC_mean.1 + 2.009*sqrt(res_df8$ANWC_var.1) < 0.3)/1000
1-sum(res_df8$ANWC_mean.2 - 2.009*sqrt(res_df8$ANWC_var.2) > -2|res_df8$ANWC_mean.2 + 2.009*sqrt(res_df8$ANWC_var.2) < -2)/1000
1-sum(res_df8$ANWC_mean.3 - 2.009*sqrt(res_df8$ANWC_var.3) > -0.25|res_df8$ANWC_mean.3 + 2.009*sqrt(res_df8$ANWC_var.3) < -0.25)/1000
1-sum(res_df8$ANWC_mean.4 - 2.009*sqrt(res_df8$ANWC_var.4) > 0.1|res_df8$ANWC_mean.4 + 2.009*sqrt(res_df8$ANWC_var.4) < 0.1)/1000
1-sum(res_df8$ANWC_mean.5 - 2.009*sqrt(res_df8$ANWC_var.5) > -1.5|res_df8$ANWC_mean.5 + 2.009*sqrt(res_df8$ANWC_var.5) < -1.5)/1000
1-sum(res_df8$ANWC_T - 2.009*sqrt(res_df8$ANWC_TV) > 48030 | res_df8$ANWC_T + 2.009*sqrt(res_df8$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df8$ANWC_var.5)
var(res_df8$ANWC_mean.5)
mean(res_df8$ANWC_var.4)
var(res_df8$ANWC_mean.4)
mean(res_df8$ANWC_var.3)
var(res_df8$ANWC_mean.3)
mean(res_df8$ANWC_var.2)
var(res_df8$ANWC_mean.2)
mean(res_df8$ANWC_var.1)
var(res_df8$ANWC_mean.1)
mean(res_df8$ANWC_TV)
var(res_df8$ANWC_T)

p <- qplot(res_df8$ANWC_mean.1,alpha_premiss7[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df8$ANWC_mean.2,alpha_premiss7[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100

save(res_df7,alpha_premiss7,TPR,TNR,FPR,FNR,file="./results7/res_df7.RData")


TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results4/Mis_",i,".RData",sep=""))
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# 03/19 Update
# process results9
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss9 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results9/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss9[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# coverage of CI
res_df9 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df9$ANWC_mean.1 - 2.009*sqrt(res_df9$ANWC_var.1) > 0.3|res_df9$ANWC_mean.1 + 2.009*sqrt(res_df9$ANWC_var.1) < 0.3)/1000
1-sum(res_df9$ANWC_mean.2 - 2.009*sqrt(res_df9$ANWC_var.2) > -0.5|res_df9$ANWC_mean.2 + 2.009*sqrt(res_df9$ANWC_var.2) < -2)/1000
1-sum(res_df9$ANWC_mean.3 - 2.009*sqrt(res_df9$ANWC_var.3) > -0.25|res_df9$ANWC_mean.3 + 2.009*sqrt(res_df9$ANWC_var.3) < -0.25)/1000
1-sum(res_df9$ANWC_mean.4 - 2.009*sqrt(res_df9$ANWC_var.4) > 0.5|res_df9$ANWC_mean.4 + 2.009*sqrt(res_df9$ANWC_var.4) < 0.1)/1000
1-sum(res_df9$ANWC_mean.5 - 2.009*sqrt(res_df9$ANWC_var.5) > 0|res_df9$ANWC_mean.5 + 2.009*sqrt(res_df9$ANWC_var.5) < -1.5)/1000
1-sum(res_df9$ANWC_T - 2.009*sqrt(res_df9$ANWC_TV) > 64095 | res_df9$ANWC_T + 2.009*sqrt(res_df9$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df9$ANWC_var.5)
var(res_df9$ANWC_mean.5)
mean(res_df9$ANWC_var.4)
var(res_df9$ANWC_mean.4)
mean(res_df9$ANWC_var.3)
var(res_df9$ANWC_mean.3)
mean(res_df9$ANWC_var.2)
var(res_df9$ANWC_mean.2)
mean(res_df9$ANWC_var.1)
var(res_df9$ANWC_mean.1)
mean(res_df9$ANWC_TV)
var(res_df9$ANWC_T)

p <- qplot(res_df9$ANWC_mean.1,alpha_premiss9[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df9$ANWC_mean.2,alpha_premiss9[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100

save(res_df9,alpha_premiss9,TPR,TNR,FPR,FNR,file="./results9/res_df9.RData")


# process results10
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss10 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results10/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss10[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
}

# coverage of CI
res_df10 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df10$ANWC_mean.1 - 2.009*sqrt(res_df10$ANWC_var.1) > 0.3|res_df10$ANWC_mean.1 + 2.009*sqrt(res_df10$ANWC_var.1) < 0.3)/1000
1-sum(res_df10$ANWC_mean.2 - 2.009*sqrt(res_df10$ANWC_var.2) > -0.5|res_df10$ANWC_mean.2 + 2.009*sqrt(res_df10$ANWC_var.2) < -2)/1000
1-sum(res_df10$ANWC_mean.3 - 2.009*sqrt(res_df10$ANWC_var.3) > -0.25|res_df10$ANWC_mean.3 + 2.009*sqrt(res_df10$ANWC_var.3) < -0.25)/1000
1-sum(res_df10$ANWC_mean.4 - 2.009*sqrt(res_df10$ANWC_var.4) > 0.5|res_df10$ANWC_mean.4 + 2.009*sqrt(res_df10$ANWC_var.4) < 0.1)/1000
1-sum(res_df10$ANWC_mean.5 - 2.009*sqrt(res_df10$ANWC_var.5) > 0|res_df10$ANWC_mean.5 + 2.009*sqrt(res_df10$ANWC_var.5) < -1.5)/1000
1-sum(res_df10$ANWC_T - 2.009*sqrt(res_df10$ANWC_TV) > 64095 | res_df10$ANWC_T + 2.009*sqrt(res_df10$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df10$ANWC_var.5)
var(res_df10$ANWC_mean.5)
mean(res_df10$ANWC_var.4)
var(res_df10$ANWC_mean.4)
mean(res_df10$ANWC_var.3)
var(res_df10$ANWC_mean.3)
mean(res_df10$ANWC_var.2)
var(res_df10$ANWC_mean.2)
mean(res_df10$ANWC_var.1)
var(res_df10$ANWC_mean.1)
mean(res_df10$ANWC_TV)
var(res_df10$ANWC_T)

p <- qplot(res_df10$ANWC_mean.1,alpha_premiss10[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df10$ANWC_mean.2,alpha_premiss10[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100

save(res_df10,alpha_premiss10,TPR,TNR,FPR,FNR,file="./results10/res_df10.RData")



# 03/20 Update: get the pre-missing estimators
# process results2:
alpha_premiss2 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
total_premiss2 <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results2/Mis_",i,".RData",sep=""))
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss2[i,] <- coef(m1)
  total_premiss2[i] <- sum(sub_dat$X1*sub_dat$W)
}
colMeans(alpha_premiss2)
mean(total_premiss2)
apply(alpha_premiss2,2,var)
var(total_premiss2)

# process results3:
alpha_premiss3 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
total_premiss3 <- rep(NA,sim_n)
for (i in 1:sim_n){
  load(paste("~/Missing/results3/Mis_",i,".RData",sep=""))
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss3[i,] <- coef(m1)
  total_premiss3[i] <- sum(sub_dat$X1*sub_dat$W)
}
colMeans(alpha_premiss3)
mean(total_premiss3)
apply(alpha_premiss3,2,var)
var(total_premiss3)

# 03/29 Update
# process results11
# process results6
sim_n <- 1000
MAR_mean <- MAR_var <- ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListMAR$alpha_mean) + length(resultListMAR$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss11 <- matrix(NA,sim_n,length(resultListMAR$alpha_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- rep(NA,sim_n)
V_adj <- 10
for (i in 1:sim_n){
  load(paste("~/Missing/results11/Mis_",i,".RData",sep=""))
  MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss11[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
}

# coverage of CI
res_df11 <- data.frame(MAR_mean = MAR_mean,MAR_var = MAR_var, ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(res_df11$ANWC_mean.1 - 2.009*sqrt(res_df11$ANWC_var.1) > 0.3|res_df11$ANWC_mean.1 + 2.009*sqrt(res_df11$ANWC_var.1) < 0.3)/1000
1-sum(res_df11$ANWC_mean.2 - 2.009*sqrt(res_df11$ANWC_var.2) > -0.5|res_df11$ANWC_mean.2 + 2.009*sqrt(res_df11$ANWC_var.2) < -1)/1000
1-sum(res_df11$ANWC_mean.3 - 2.009*sqrt(res_df11$ANWC_var.3) > -0.25|res_df11$ANWC_mean.3 + 2.009*sqrt(res_df11$ANWC_var.3) < -0.25)/1000
1-sum(res_df11$ANWC_mean.4 - 2.009*sqrt(res_df11$ANWC_var.4) > 0.1|res_df11$ANWC_mean.4 + 2.009*sqrt(res_df11$ANWC_var.4) < 0.1)/1000
1-sum(res_df11$ANWC_mean.5 - 2.009*sqrt(res_df11$ANWC_var.5) > -1.5|res_df11$ANWC_mean.5 + 2.009*sqrt(res_df11$ANWC_var.5) < -1.5)/1000
1-sum(res_df11$ANWC_T - 2.009*sqrt(res_df11$ANWC_TV) > 64110 | res_df11$ANWC_T + 2.009*sqrt(res_df11$ANWC_TV) < 48030)/1000

# MCMC variance
mean(res_df11$ANWC_var.5)
var(res_df11$ANWC_mean.5)
mean(res_df11$ANWC_var.4)
var(res_df11$ANWC_mean.4)
mean(res_df11$ANWC_var.3)
var(res_df11$ANWC_mean.3)
mean(res_df11$ANWC_var.2)
var(res_df11$ANWC_mean.2)
mean(res_df11$ANWC_var.1)
var(res_df11$ANWC_mean.1)
mean(res_df11$ANWC_TV)
var(res_df11$ANWC_T)

p <- qplot(res_df11$ANWC_mean.1,alpha_premiss11[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(res_df11$ANWC_mean.2,alpha_premiss11[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100

