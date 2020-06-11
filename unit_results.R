# 04/20 Update
# process unit1

sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + 
                                  length(resultListANWC$gamma_mean)+length(resultListANWC$beta_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_premiss1 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 10
for (i in 1:sim_n){
  load(paste("~/Missing/unit1/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
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
  alpha_premiss1[i,] <- coef(m1)
  mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
}

# coverage of CI
unit_df1 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var)
1-sum(unit_df1$ANWC_mean.1 - 2.009*sqrt(unit_df1$ANWC_var.1) > 0.3|unit_df1$ANWC_mean.1 + 2.009*sqrt(unit_df1$ANWC_var.1) < 0.3)/1000
1-sum(unit_df1$ANWC_mean.2 - 2.009*sqrt(unit_df1$ANWC_var.2) > -0.5|unit_df1$ANWC_mean.2 + 2.009*sqrt(unit_df1$ANWC_var.2) < -0.5)/1000
1-sum(unit_df1$ANWC_mean.3 - 2.009*sqrt(unit_df1$ANWC_var.3) > -0.25|unit_df1$ANWC_mean.3 + 2.009*sqrt(unit_df1$ANWC_var.3) < -0.25)/1000
1-sum(unit_df1$ANWC_mean.4 - 2.009*sqrt(unit_df1$ANWC_var.4) > 0.1|unit_df1$ANWC_mean.4 + 2.009*sqrt(unit_df1$ANWC_var.4) < 0.1)/1000
1-sum(unit_df1$ANWC_mean.5 - 2.009*sqrt(unit_df1$ANWC_var.5) > -1.5|unit_df1$ANWC_mean.5 + 2.009*sqrt(unit_df1$ANWC_var.5) < -1.5)/1000
1-sum(unit_df1$ANWC_T - 2.009*sqrt(unit_df1$ANWC_TV) >  64156 | unit_df1$ANWC_T + 2.009*sqrt(unit_df1$ANWC_TV) <  64156)/1000

# MCMC variance
mean(unit_df1$ANWC_var.5)
var(unit_df1$ANWC_mean.5)
mean(unit_df1$ANWC_var.4)
var(unit_df1$ANWC_mean.4)
mean(unit_df1$ANWC_var.3)
var(unit_df1$ANWC_mean.3)
mean(unit_df1$ANWC_var.2)
var(unit_df1$ANWC_mean.2)
mean(unit_df1$ANWC_var.1)
var(unit_df1$ANWC_mean.1)
mean(unit_df1$ANWC_TV)
var(unit_df1$ANWC_T)

p <- qplot(unit_df1$ANWC_mean.1,alpha_premiss11[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(unit_df1$ANWC_mean.2,alpha_premiss11[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100


# 05/08 update
# process unit2
sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_allunit2 <- alpha_premiss2 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- Ty_total_var <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 1
for (i in 1:sim_n){
  load(paste("~/Missing/unit2/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  mydesign <- svydesign(id = ~id,data = test_dat[which(test_dat$U==0),],weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_allunit2[i,] <- coef(m1)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss2[i,] <- coef(m2)
  #mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  #imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  #TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  #TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  #FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  #FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
  Ty_total_var[i] <- sum((as.numeric(as.character(test_dat$Y))*test_dat$W)^2*(1-1/test_dat$W))
}

# coverage of CI
unit_df2 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var,ANWC_Y = Ty,ANWC_Y_var = Ty_total_var)
1-sum(unit_df2$ANWC_mean.1 - 2.009*sqrt(unit_df2$ANWC_var.1) > 0.3|unit_df2$ANWC_mean.1 + 2.009*sqrt(unit_df2$ANWC_var.1) < 0.3)/1000
1-sum(unit_df2$ANWC_mean.2 - 2.009*sqrt(unit_df2$ANWC_var.2) > -0.5|unit_df2$ANWC_mean.2 + 2.009*sqrt(unit_df2$ANWC_var.2) < -0.5)/1000
1-sum(unit_df2$ANWC_mean.3 - 2.009*sqrt(unit_df2$ANWC_var.3) > -0.25|unit_df2$ANWC_mean.3 + 2.009*sqrt(unit_df2$ANWC_var.3) < -0.25)/1000
1-sum(unit_df2$ANWC_mean.4 - 2.009*sqrt(unit_df2$ANWC_var.4) > 0.1|unit_df2$ANWC_mean.4 + 2.009*sqrt(unit_df2$ANWC_var.4) < 0.1)/1000
1-sum(unit_df2$ANWC_mean.5 - 2.009*sqrt(unit_df2$ANWC_var.5) > -1.5|unit_df2$ANWC_mean.5 + 2.009*sqrt(unit_df2$ANWC_var.5) < -1.5)/1000
1-sum(unit_df2$ANWC_T - 2.009*sqrt(unit_df2$ANWC_TV) >  64156 | unit_df2$ANWC_T + 2.009*sqrt(unit_df2$ANWC_TV) <  64156)/1000
1-sum(unit_df2$ANWC_Y - 2.009*sqrt(unit_df2$ANWC_Y_var) >  sum(pop_dat$Y) | unit_df2$ANWC_Y + 2.009*sqrt(unit_df2$ANWC_Y_var) <  sum(pop_dat$Y))/1000

# MCMC variance
mean(unit_df2$ANWC_var.5)
var(unit_df2$ANWC_mean.5)
mean(unit_df2$ANWC_var.4)
var(unit_df2$ANWC_mean.4)
mean(unit_df2$ANWC_var.3)
var(unit_df2$ANWC_mean.3)
mean(unit_df2$ANWC_var.2)
var(unit_df2$ANWC_mean.2)
mean(unit_df2$ANWC_var.1)
var(unit_df2$ANWC_mean.1)
mean(unit_df2$ANWC_TV)
var(unit_df2$ANWC_T)
mean(unit_df2$ANWC_Y_var)
var(unit_df2$ANWC_Y)

p <- qplot(unit_df2$ANWC_mean.1,alpha_premiss2[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)
p <- qplot(unit_df2$ANWC_mean.1,alpha_allunit2[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "all unit")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(unit_df2$ANWC_mean.2,alpha_premiss2[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)
p <- qplot(unit_df2$ANWC_mean.2,alpha_allunit2[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "all unit")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100


# 05/09 update
# process unit3
sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_allunit3 <- alpha_premiss3 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- Ty_total_var <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 1
for (i in 1:sim_n){
  load(paste("~/Missing/unit3/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  mydesign <- svydesign(id = ~id,data = test_dat[which(test_dat$U==0),],weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_allunit3[i,] <- coef(m1)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss3[i,] <- coef(m2)
  #mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  #imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  #TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  #TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  #FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  #FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
  Ty_total_var[i] <- sum((as.numeric(as.character(test_dat$Y))*test_dat$W)^2*(1-1/test_dat$W))
}

# coverage of CI
unit_df3 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var,ANWC_Y = Ty,ANWC_Y_var = Ty_total_var)
1-sum(unit_df3$ANWC_mean.1 - 2.009*sqrt(unit_df3$ANWC_var.1) > 0.3|unit_df3$ANWC_mean.1 + 2.009*sqrt(unit_df3$ANWC_var.1) < 0.3)/1000
1-sum(unit_df3$ANWC_mean.2 - 2.009*sqrt(unit_df3$ANWC_var.2) > -0.5|unit_df3$ANWC_mean.2 + 2.009*sqrt(unit_df3$ANWC_var.2) < -0.5)/1000
1-sum(unit_df3$ANWC_mean.3 - 2.009*sqrt(unit_df3$ANWC_var.3) > -0.25|unit_df3$ANWC_mean.3 + 2.009*sqrt(unit_df3$ANWC_var.3) < -0.25)/1000
1-sum(unit_df3$ANWC_mean.4 - 2.009*sqrt(unit_df3$ANWC_var.4) > 0.1|unit_df3$ANWC_mean.4 + 2.009*sqrt(unit_df3$ANWC_var.4) < 0.1)/1000
1-sum(unit_df3$ANWC_mean.5 - 2.009*sqrt(unit_df3$ANWC_var.5) > -1.5|unit_df3$ANWC_mean.5 + 2.009*sqrt(unit_df3$ANWC_var.5) < -1.5)/1000
1-sum(unit_df3$ANWC_T - 2.009*sqrt(unit_df3$ANWC_TV) >  64156 | unit_df3$ANWC_T + 2.009*sqrt(unit_df3$ANWC_TV) <  64156)/1000
1-sum(unit_df3$ANWC_Y - 2.009*sqrt(unit_df3$ANWC_Y_var) >  sum(pop_dat$Y) | unit_df3$ANWC_Y + 2.009*sqrt(unit_df3$ANWC_Y_var) <  sum(pop_dat$Y))/1000

# MCMC variance
mean(unit_df3$ANWC_var.5)
var(unit_df3$ANWC_mean.5)
mean(unit_df3$ANWC_var.4)
var(unit_df3$ANWC_mean.4)
mean(unit_df3$ANWC_var.3)
var(unit_df3$ANWC_mean.3)
mean(unit_df3$ANWC_var.2)
var(unit_df3$ANWC_mean.2)
mean(unit_df3$ANWC_var.1)
var(unit_df3$ANWC_mean.1)
mean(unit_df3$ANWC_TV)
var(unit_df3$ANWC_T)
mean(unit_df3$ANWC_Y_var)
var(unit_df3$ANWC_Y)

p <- qplot(unit_df3$ANWC_mean.1,alpha_premiss3[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)
p <- qplot(unit_df3$ANWC_mean.1,alpha_allunit3[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "all unit")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(unit_df3$ANWC_mean.2,alpha_premiss3[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)
p <- qplot(unit_df3$ANWC_mean.2,alpha_allunit3[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "all unit")
p + geom_abline(intercept = 0, slope = 1)

mean(TPR)*100
mean(TNR)*100
mean(FPR)*100
mean(FNR)*100
sd(TPR)/sqrt(1000)*100
sd(TNR)/sqrt(1000)*100
sd(FPR)/sqrt(1000)*100
sd(FNR)/sqrt(1000)*100

save(unit_df3,alpha_allunit3,alpha_premiss3,file = "./unit3/Results.RData")

# 05/10 update
# process unit4
sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_allunit4 <- alpha_premiss4 <- beta_allunit4 <- beta_premiss4 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- Ty_total_var <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 1
for (i in 1:sim_n){
  load(paste("~/Missing/unit4/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  mydesign <- svydesign(id = ~id,data = test_dat[which(test_dat$U==0),],weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  m3 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_allunit4[i,] <- coef(m1)
  beta_allunit4[i,] <- coef(m3)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss4[i,] <- coef(m2)
  beta_premiss4[i,] <- coef(m4)
  #mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  #imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  #TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  #TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  #FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  #FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
  Ty_total_var[i] <- sum((as.numeric(as.character(test_dat$Y))*test_dat$W)^2*(1-1/test_dat$W))
}

# coverage of CI
unit_df4 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var,ANWC_Y = Ty,ANWC_Y_var = Ty_total_var)
1-sum(unit_df4$ANWC_mean.1 - 2.009*sqrt(unit_df4$ANWC_var.1) > 0.3|unit_df4$ANWC_mean.1 + 2.009*sqrt(unit_df4$ANWC_var.1) < 0.3)/1000
1-sum(unit_df4$ANWC_mean.2 - 2.009*sqrt(unit_df4$ANWC_var.2) > -0.5|unit_df4$ANWC_mean.2 + 2.009*sqrt(unit_df4$ANWC_var.2) < -0.5)/1000
1-sum(unit_df4$ANWC_mean.3 - 2.009*sqrt(unit_df4$ANWC_var.3) > -0.25|unit_df4$ANWC_mean.3 + 2.009*sqrt(unit_df4$ANWC_var.3) < -0.25)/1000
1-sum(unit_df4$ANWC_mean.4 - 2.009*sqrt(unit_df4$ANWC_var.4) > 0.1|unit_df4$ANWC_mean.4 + 2.009*sqrt(unit_df4$ANWC_var.4) < 0.1)/1000
1-sum(unit_df4$ANWC_mean.5 - 2.009*sqrt(unit_df4$ANWC_var.5) > -1.5|unit_df4$ANWC_mean.5 + 2.009*sqrt(unit_df4$ANWC_var.5) < -1.5)/1000
1-sum(unit_df4$ANWC_T - 2.009*sqrt(unit_df4$ANWC_TV) >  64156 | unit_df4$ANWC_T + 2.009*sqrt(unit_df4$ANWC_TV) <  64156)/1000
1-sum(unit_df4$ANWC_Y - 2.009*sqrt(unit_df4$ANWC_Y_var) >  sum(pop_dat$Y) | unit_df4$ANWC_Y + 2.009*sqrt(unit_df4$ANWC_Y_var) <  sum(pop_dat$Y))/1000

# MCMC variance
mean(unit_df4$ANWC_var.5)
var(unit_df4$ANWC_mean.5)
mean(unit_df4$ANWC_var.4)
var(unit_df4$ANWC_mean.4)
mean(unit_df4$ANWC_var.3)
var(unit_df4$ANWC_mean.3)
mean(unit_df4$ANWC_var.2)
var(unit_df4$ANWC_mean.2)
mean(unit_df4$ANWC_var.1)
var(unit_df4$ANWC_mean.1)
mean(unit_df4$ANWC_TV)
var(unit_df4$ANWC_T)
mean(unit_df4$ANWC_Y_var)
var(unit_df4$ANWC_Y)

p <- qplot(unit_df4$ANWC_mean.1,alpha_premiss4[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)
p <- qplot(unit_df4$ANWC_mean.1,alpha_allunit4[,1],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "all unit")
p + geom_abline(intercept = 0, slope = 1)

p <- qplot(unit_df4$ANWC_mean.2,alpha_premiss4[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "pre-missing")
p + geom_abline(intercept = 0, slope = 1)
p <- qplot(unit_df4$ANWC_mean.2,alpha_allunit4[,2],geom="point", alpha=I(0.5),xlab = "estimate",ylab = "all unit")
p + geom_abline(intercept = 0, slope = 1)


save(unit_df4,alpha_allunit4,alpha_premiss4,file = "./unit4/Results.RData")


# 05/14 update
# process unit4
sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + length(resultListANWC$gamma_mean))
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_allunit5 <- alpha_premiss5 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- Ty_total_var <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 1
for (i in 1:sim_n){
  load(paste("~/Missing/unit5/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  mydesign <- svydesign(id = ~id,data = test_dat[which(test_dat$U==0),],weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_allunit5[i,] <- coef(m1)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss5[i,] <- coef(m2)
  #mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  #imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  #TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  #TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  #FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  #FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
  Ty_total_var[i] <- sum((as.numeric(as.character(test_dat$Y))*test_dat$W)^2*(1-1/test_dat$W))
}

# point estimates 
mean(unit_df5$ANWC_mean.1)
mean(unit_df5$ANWC_mean.2)
mean(unit_df5$ANWC_mean.3)
mean(unit_df5$ANWC_mean.4)
mean(unit_df5$ANWC_mean.5)

# examine imputed Y
library(dplyr)
for (i in 1:sim_n){
  load(paste("~/Missing/unit5/Mis_",i,".RData",sep=""))
  trueY <- sub_dat$Y
  imputeY <- t(select(MI_dataANWC,contains("Y")))
  imputeX <- t(select(MI_dataANWC,contains("X")))
}


# 05/24 update
# process unit7
sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + 
                                  length(resultListANWC$gamma_mean)+2)
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_allunit7 <- alpha_premiss7 <- beta_allunit7 <- beta_premiss7 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
gamma_allunit7 <- gamma_premiss7 <- matrix(NA,sim_n,length(resultListANWC$gamma_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- Ty_total_var <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 1
for (i in 1:sim_n){
  load(paste("~/Missing/unit7/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  
  # colMeans(resultListANWC$beta) should be resultListANWC$beta after modification of code
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean,resultListANWC$beta_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var,resultListANWC$beta_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat[which(test_dat$U==0),],weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  m3 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
  m5 <- svyglm(Rx~Y+X1,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_allunit7[i,] <- coef(m1)
  beta_allunit7[i,] <- coef(m3)
  gamma_allunit7[i,] <- coef(m5)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
  m6 <- svyglm(Rx~Y+X1,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss7[i,] <- coef(m2)
  beta_premiss7[i,] <- coef(m4)
  gamma_premiss7[i,] <- coef(m6)
  #mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  #imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  #TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  #TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  #FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  #FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
  Ty_total_var[i] <- sum((as.numeric(as.character(test_dat$Y))*test_dat$W)^2*(1-1/test_dat$W))
}

# coverage of CI
unit_df7 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var,ANWC_Y = Ty,ANWC_Y_var = Ty_total_var)
1-sum(unit_df7$ANWC_mean.1 - 2.009*sqrt(unit_df7$ANWC_var.1) > 0.3|unit_df7$ANWC_mean.1 + 2.009*sqrt(unit_df7$ANWC_var.1) < 0.3)/1000
1-sum(unit_df7$ANWC_mean.2 - 2.009*sqrt(unit_df7$ANWC_var.2) > -0.5|unit_df7$ANWC_mean.2 + 2.009*sqrt(unit_df7$ANWC_var.2) < -0.5)/1000
1-sum(unit_df7$ANWC_mean.3 - 2.009*sqrt(unit_df7$ANWC_var.3) > -0.25|unit_df7$ANWC_mean.3 + 2.009*sqrt(unit_df7$ANWC_var.3) < -0.25)/1000
1-sum(unit_df7$ANWC_mean.4 - 2.009*sqrt(unit_df7$ANWC_var.4) > 0.1|unit_df7$ANWC_mean.4 + 2.009*sqrt(unit_df7$ANWC_var.4) < 0.1)/1000
1-sum(unit_df7$ANWC_mean.5 - 2.009*sqrt(unit_df7$ANWC_var.5) > -1.5|unit_df7$ANWC_mean.5 + 2.009*sqrt(unit_df7$ANWC_var.5) < -1.5)/1000
1-sum(unit_df7$ANWC_mean.6 - 2.009*sqrt(unit_df7$ANWC_var.6) > 0.5|unit_df7$ANWC_mean.6 + 2.009*sqrt(unit_df7$ANWC_var.6) < 0.5)/1000
1-sum(unit_df7$ANWC_mean.7 - 2.009*sqrt(unit_df7$ANWC_var.7) > 0.015|unit_df7$ANWC_mean.7 + 2.009*sqrt(unit_df7$ANWC_var.7) < 0.015)/1000
1-sum(unit_df7$ANWC_T - 2.009*sqrt(unit_df7$ANWC_TV) >  64156 | unit_df7$ANWC_T + 2.009*sqrt(unit_df7$ANWC_TV) <  64156)/1000
1-sum(unit_df7$ANWC_Y - 2.009*sqrt(unit_df7$ANWC_Y_var) >  sum(pop_dat$Y) | unit_df7$ANWC_Y + 2.009*sqrt(unit_df7$ANWC_Y_var) <  sum(pop_dat$Y))/1000

# MCMC variance
mean(unit_df7$ANWC_var.7)
var(unit_df7$ANWC_mean.7)
mean(unit_df7$ANWC_var.6)
var(unit_df7$ANWC_mean.6)
mean(unit_df7$ANWC_var.5)
var(unit_df7$ANWC_mean.5)
mean(unit_df7$ANWC_var.4)
var(unit_df7$ANWC_mean.4)
mean(unit_df7$ANWC_var.3)
var(unit_df7$ANWC_mean.3)
mean(unit_df7$ANWC_var.2)
var(unit_df7$ANWC_mean.2)
mean(unit_df7$ANWC_var.1)
var(unit_df7$ANWC_mean.1)
mean(unit_df7$ANWC_TV)
var(unit_df7$ANWC_T)
mean(unit_df7$ANWC_Y_var)
var(unit_df7$ANWC_Y)

# 05/26 update
# process unit8
sim_n <- 1000
ANWC_mean <- ANWC_var <- matrix(NA,sim_n,length(resultListANWC$alpha_mean) + 
                                  length(resultListANWC$gamma_mean)+2)
ANWC_total <- ANWC_total_var <- rep(NA,sim_n)
alpha_allunit8 <- alpha_premiss8 <- beta_allunit8 <- beta_premiss8 <- matrix(NA,sim_n,length(resultListANWC$alpha_mean))
gamma_allunit8 <- gamma_premiss8 <- matrix(NA,sim_n,length(resultListANWC$gamma_mean))
TPR <- TNR <- FPR <- FNR <- pop_sd_MAR <- pop_sd <- Ty <- Ty_total_var <- rep(NA,sim_n)
accpet_ratio <- rep(NA,sim_n)
V_adj <- 1
for (i in 1:sim_n){
  load(paste("~/Missing/unit8/Mis_",i,".RData",sep=""))
  #MAR_mean[i,] <- c(resultListMAR$alpha_mean, resultListMAR$gamma_mean)
  #MAR_var[i,] <- c(resultListMAR$alpha_var, resultListMAR$gamma_var)
  
  # colMeans(resultListANWC$beta) should be resultListANWC$beta after modification of code
  ANWC_mean[i,] <- c(resultListANWC$alpha_mean, resultListANWC$gamma_mean,resultListANWC$beta_mean)
  ANWC_var[i,] <- c(resultListANWC$alpha_var, resultListANWC$gamma_var,resultListANWC$beta_var)
  ANWC_total[i] <- resultListANWC$total_mean
  ANWC_total_var[i] <- resultListANWC$total_var
  test_dat <- as.data.frame(cbind(sub_dat$X1,rownames(sub_dat),sub_dat$W,sub_dat$Y,sub_dat$U,sub_dat$Rx),stringsAsFactors = FALSE)
  names(test_dat) <- c("X1","id","W","Y","U","Rx")
  test_dat$X1 <- as.factor(test_dat$X1)
  test_dat$id <- as.numeric(test_dat$id)
  test_dat$W <- as.numeric(test_dat$W)
  test_dat$pi <- 1/test_dat$W
  test_dat$Y <- as.factor(test_dat$Y)
  test_dat$U <- as.factor(test_dat$U)
  test_dat$Rx <- as.factor(test_dat$Rx)
  mydesign <- svydesign(id = ~id,data = test_dat[which(test_dat$U==0),],weight = ~W)
  m1 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  m3 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
  m5 <- svyglm(Rx~Y+X1,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_allunit8[i,] <- coef(m1)
  beta_allunit8[i,] <- coef(m3)
  gamma_allunit8[i,] <- coef(m5)
  mydesign <- svydesign(id = ~id,data = test_dat,weight = ~W)
  m2 <- svyglm(X1~Y,design = mydesign,family = quasibinomial(link = "probit"))
  m4 <- svyglm(Y~pi,design = mydesign,family = quasibinomial(link = "probit"))
  m6 <- svyglm(Rx~Y+X1,design = mydesign,family = quasibinomial(link = "probit"))
  alpha_premiss8[i,] <- coef(m2)
  beta_premiss8[i,] <- coef(m4)
  gamma_premiss8[i,] <- coef(m6)
  #mis_datX <- sub_dat[which(sub_dat$Rx == 1),]$X1
  #imp_dat <- MI_dataANWC[,which(sub_dat$Rx == 1)]
  #TPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 1)/length(mis_datX) # truly 1, imputed 1
  #TNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 0)/length(mis_datX)
  #FPR[i] <- sum(imp_dat[1,] == 1 & mis_datX == 0)/length(mis_datX) # imputed 1, but 0 in truth
  #FNR[i] <- sum(imp_dat[1,] == 0 & mis_datX == 1)/length(mis_datX)
  pop_sd_MAR[i] <- pop_sd_HT
  pop_sd[i] <- sqrt(sum((sub_dat$X1/sub_dat$pi)^2*(1-sub_dat$pi)))/V_adj
  accpet_ratio[i] <- ratio
  Ty[i] <- sum(as.numeric(as.character((test_dat$Y)))*test_dat$W)
  Ty_total_var[i] <- sum((as.numeric(as.character(test_dat$Y))*test_dat$W)^2*(1-1/test_dat$W))
}

# coverage of CI
unit_df8 <- data.frame(ANWC_mean = ANWC_mean, ANWC_var = ANWC_var, ANWC_T = ANWC_total, ANWC_TV = ANWC_total_var,ANWC_Y = Ty,ANWC_Y_var = Ty_total_var)
1-sum(unit_df8$ANWC_mean.1 - 2.009*sqrt(unit_df8$ANWC_var.1) > 0.3|unit_df8$ANWC_mean.1 + 2.009*sqrt(unit_df8$ANWC_var.1) < 0.3)/1000
1-sum(unit_df8$ANWC_mean.2 - 2.009*sqrt(unit_df8$ANWC_var.2) > -0.5|unit_df8$ANWC_mean.2 + 2.009*sqrt(unit_df8$ANWC_var.2) < -0.5)/1000
1-sum(unit_df8$ANWC_mean.3 - 2.009*sqrt(unit_df8$ANWC_var.3) > -0.25|unit_df8$ANWC_mean.3 + 2.009*sqrt(unit_df8$ANWC_var.3) < -0.25)/1000
1-sum(unit_df8$ANWC_mean.4 - 2.009*sqrt(unit_df8$ANWC_var.4) > 0.1|unit_df8$ANWC_mean.4 + 2.009*sqrt(unit_df8$ANWC_var.4) < 0.1)/1000
1-sum(unit_df8$ANWC_mean.5 - 2.009*sqrt(unit_df8$ANWC_var.5) > -1.5|unit_df8$ANWC_mean.5 + 2.009*sqrt(unit_df8$ANWC_var.5) < -1.5)/1000
1-sum(unit_df8$ANWC_mean.6 - 2.009*sqrt(unit_df8$ANWC_var.6) > 0.5|unit_df8$ANWC_mean.6 + 2.009*sqrt(unit_df8$ANWC_var.6) < 0.5)/1000
1-sum(unit_df8$ANWC_mean.7 - 2.009*sqrt(unit_df8$ANWC_var.7) > 0.015|unit_df8$ANWC_mean.7 + 2.009*sqrt(unit_df8$ANWC_var.7) < 0.015)/1000
1-sum(unit_df8$ANWC_T - 2.009*sqrt(unit_df8$ANWC_TV) >  64156 | unit_df8$ANWC_T + 2.009*sqrt(unit_df8$ANWC_TV) <  64156)/1000
1-sum(unit_df8$ANWC_Y - 2.009*sqrt(unit_df8$ANWC_Y_var) >  sum(pop_dat$Y) | unit_df8$ANWC_Y + 2.009*sqrt(unit_df8$ANWC_Y_var) <  sum(pop_dat$Y))/1000

# MCMC variance
mean(unit_df8$ANWC_var.7)
var(unit_df8$ANWC_mean.7)
mean(unit_df8$ANWC_var.6)
var(unit_df8$ANWC_mean.6)
mean(unit_df8$ANWC_var.5)
var(unit_df8$ANWC_mean.5)
mean(unit_df8$ANWC_var.4)
var(unit_df8$ANWC_mean.4)
mean(unit_df8$ANWC_var.3)
var(unit_df8$ANWC_mean.3)
mean(unit_df8$ANWC_var.2)
var(unit_df8$ANWC_mean.2)
mean(unit_df8$ANWC_var.1)
var(unit_df8$ANWC_mean.1)
mean(unit_df8$ANWC_TV)
var(unit_df8$ANWC_T)
mean(unit_df8$ANWC_Y_var)
var(unit_df8$ANWC_Y)
