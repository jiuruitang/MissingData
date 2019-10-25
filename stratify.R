# Implement Michael's stratified sampling

N = 50000
N1 <- N*0.7
N2 <- N*0.3
n <- 5000
n1 <- n*0.3
n2 <- n*0.7

# 1 to 35000 people in stratum 1, 35001 to 50000 people in stratum 2
# theta = (0.5,0.15,0.35) stratum 1, theta = (0.1,0.45,0.45) stratum 2

# sample Y for multinomial
Y <- rep(NA,N)
set.seed(1233)
Y[1:N1] <- sample(c(1,2,3),N1,replace = TRUE, prob = c(0.5,0.15,0.35))
set.seed(2324)
Y[(N1+1):N] <- sample(c(1,2,3),N2,replace = TRUE, prob = c(0.1,0.45,0.45))
W <- c(rep(N1/n1,N1),rep(N2/n2,N2))

# scenario 1
alpha0 <- 0.5
alpha12 <- -0.5
alpha13 <- -1
gamma0 <- -0.25
gamma12 <- 0.1
gamma13 <- 0.3
gamma2 <- -1.1

# sample X1|Y
set.seed(45343)
pi_x1 <- pnorm(alpha0 + alpha12*(Y == 2) + alpha13*(Y==3)) # pnorm is CDF or normal
set.seed(3204)
X1 <- rbinom(N,1,p=pi_x1)

# sample Rx|X,Y
set.seed(894)
pi_Rx <- pnorm(gamma0 + gamma12*(Y == 2) + gamma13*(Y==3) + gamma2*X)
set.seed(3901)
Rx <- binom(N,1,p = pi_Rx)
