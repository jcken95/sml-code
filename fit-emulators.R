## this is the one to use

library(rstan)
source("code/GPfunctionsOptim.R")
source("code/hetGPfunctions.R")
hetGP <- stan_model("code/hetGP.stan")
sml <- stan_model("code/het-SML.stan")


## load data
#y data
y.test1 <- read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/y-val.txt",sep=",")[,1]
y.hetgp1<- read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/y-hetgp.txt",sep=",")[,1]
y.c1 <- read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/y-c.txt",sep=",")[,1]
y.e1 <- read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/y-e.txt",sep=",")[,1]

#x data
x.test <- as.matrix(read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/x-val.txt",sep=","))
x.hetgp <- as.matrix(as.matrix(read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/x-hetgp.txt",sep=",")))
x.c <- as.matrix(read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/x-c.txt",sep=","))
x.e <- as.matrix(read.table(file = "~/Documents/sensitivity-analysis/9d-emulator/emulator-data/x-e.txt",sep=","))

## plots
par(mfrow=c(3,3))
for(i in 1:9) plot(x.test[,i],y.test1,pch=20,cex=0.6)
for(i in 1:9) plot(x.hetgp[,i],y.hetgp1,pch=20,cex=0.6)
for(i in 1:9) plot(x.e[,i],y.e1,pch=20,cex=0.6)
for(i in 1:9) plot(x.c[,i],y.c1,pch=20,cex=0.6)
## define probit
## because availability in (0,1)

probit <- function(p){
  qnorm(p)
}
inv.probit <- function(q){
  pnorm(q)
}

y.test <- probit(y.test1)
y.hetgp <- probit(y.hetgp1)
y.c <- probit(y.c1)
y.e <- probit(y.e1)


par(mfrow=c(3,3))
for(i in 1:9) plot(x.test[,i],y.test,pch=20,cex=0.6)
for(i in 1:9) plot(x.hetgp[,i],y.hetgp,pch=20,cex=0.6)
for(i in 1:9) plot(x.e[,i],y.e,pch=20,cex=0.6)
for(i in 1:9) plot(x.c[,i],y.c,pch=20,cex=0.6)
## define test metrics

MSE <- function(true, pred){
  mean((true-pred)^2)
}
Score <- function(y, m, v){
  -((y-m)^2)/v - log(v)
}




## standardise inputs (x-mu)/sig

##training data
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.std <- scale(x.hetgp)

##scale validation data accordingly
x.v2 <- scale(x.test, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
x.v1 <- scale(x.test, center = attr(x.std, "scaled:center"), scale = attr(x.std, "scaled:scale"))

#### fit emulators! ####

### first fit HetGP emulator ####

data.hetGP <- list(
  m_p = 10, v_p = 10,
  N = length(y.hetgp), K = 9,
  x = x.std, 
  m_H = cbind(1, log(x.hetgp)), v_H = cbind(1, x.std),
  y = as.vector(y.hetgp),
  a = rep(1,length(y.hetgp)),
  ## prior
  m_beta_m = rep(0, 10), m_beta_s = rep(10, 10),
  m_a_theta = rep(2,9), m_b_theta = rep(1,9),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget = 0,
  v_beta_m = rep(0, 10), v_beta_s = rep(10, 10),
  v_a_theta = rep(2,9), v_b_theta = rep(1,9),
  v_a_sigma = 2, v_b_sigma = 2,
  v_nugget_a = 2, v_nugget_b = 2
)
temp <- list()

find.mode <- function(x){
  rstan::optimizing(hetGP, data = data.hetGP, verbose = F, as_vector = F)
}

temp <- parallel::mclapply(1:3, find.mode, mc.cores = 3)
beepr::beep()
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))


het.fit <- temp[[best.emulator]]


pars <- het.fit$par

## predict response for hetgp
mu.y <- cbind(1, log(x.test))%*%pars$m_beta
mu.yobs <- cbind(1, log(x.hetgp))%*%pars$m_beta

mu.y.var <- as.vector(cbind(1, x.v1)%*%pars$v_beta)
mu.yobs.var <- as.vector(cbind(1, x.std)%*%pars$v_beta)

Mean.het <- cond.mean.het(y.hetgp, x.std, x.v1, pars$m_sigma, pars$m_theta, exp(pars$logLambda), 0, mu.y, mu.yobs, 2)
Var.het <- cond.var.het(pars$logLambda, x.std, x.v1, pars$v_sigma, pars$v_theta, pars$v_nugget,mu.y.var, mu.yobs.var, 2)
Var.het <- cond.covmat.het(exp(pars$logLambda), Var.het, x.std, x.v1, pars$v_sigma, pars$v_theta, 2)

#### fit SML emulator ####

## need to get the data into shape!
n.ml <- length(y.e)

y.c2 <- rev(y.c)
y.e2 <- rev(y.e)
n.c <- length(y.c2); n.e <- length(y.e2)
x.c.std <- x.c.std[n.c:1,]
x.e.std <- x.e.std[n.ml:1,]

pairs(cbind(y.c2, x.c.std))
pairs(cbind(y.e2, x.e.std))
## reverse everything to get the order correct ...
n.c <- length(y.c2); n.e <- length(y.e2)

## indexing "flips" the matrix upside down so that the design matrix is properly "aligned"

m.h <- cbind(cbind(1, log(x.c[n.c:1,])), rbind(matrix(0, ncol = 10, nrow = n.c - n.e) ,cbind(1, log(x.e[n.e:1,]))))
tail(m.h)
ml.data <- list(
  
  ## data ##
  
  m_p = 10, v_p = 10,
  N = n.e + n.c, K = 9,
  n_c = n.c, n_e = n.e,
  x_e = x.e.std, x_c = x.c.std,
  m_H = m.h, v_H = cbind(1, x.e.std),
  y_c = y.c2, y_e = y.e2,
  
  ## priors ##
  
  m_beta_m = rep(0, 10), m_beta_s = rep(10, 10),
  m_a_theta = rep(2, 9), m_b_theta = rep(1, 9),
  m_a_sigma = 2, m_b_sigma = 2,
  m_nugget = 0,
  
  v_beta_m = rep(0, 10), v_beta_s = rep(10, 10),
  v_a_theta = rep(2, 9), v_b_theta = rep(1, 9),
  v_a_sigma = 2, v_b_sigma = 2,
  v_nugget_a = 2, v_nugget_b = 2,
  
  c_beta_m = rep(0, 10), c_beta_s = rep(10, 10),
  c_a_theta = rep(2, 9), c_b_theta = rep(1, 9),
  c_a_sigma = 2, c_b_sigma = 2,
  c_nugget_a = 2, c_nugget_b = 2,
  m_rho = 1, s_rho = 1/3
  
)


temp <- list()

## fit multilevel GP

find.mode <- function(x){
  rstan::optimizing(sml, data = ml.data, verbose = F, as_vector = F)
}
st1 <- Sys.time()
temp <- parallel::mclapply(1:3, find.mode, mc.cores = 3)
en1 <- Sys.time()
en1-st1
beepr::beep(4)
best.emulator <- which.max(c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value))
c(temp[[1]]$value , temp[[2]]$value, temp[[3]]$value)

ml.fit <-  temp[[best.emulator]]
pars <- ml.fit$par


mu.yobs <- c(cbind(1,log(x.c[n.c:1,]))%*%pars$c_beta, cbind(1,log(x.e[n.e:1,]))%*%(pars$rho*pars$c_beta + pars$m_beta))
mu.y <- cbind(1,log(x.test))%*%(pars$rho*pars$c_beta + pars$m_beta)

Mean.ml <- cond.mean.hetSML(c(y.c2, y.e2), rbind(x.c.std, x.e.std), x.v2,
                            c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                            c(pars$m_sigma, pars$m_theta, 0), exp(pars$logLambda), pars$rho,
                            mu.y, mu.yobs, n.c, 2)
Var.ml <- cond.var.het(pars$logLambda, x.e.std, x.v2,
                       pars$v_sigma, pars$v_theta, pars$v_nugget,
                       as.vector(cbind(1,x.v2)%*%pars$v_beta), 
                       as.vector(cbind(1, x.e.std)%*%pars$v_beta), 2)
Var.ml <-  cond.covmat.hetSML( exp(pars$logLambda), Var.ml, rbind(x.c.std, x.e.std), x.v2,
                               pars$m_sigma, pars$m_theta, pars$c_sigma, pars$c_theta, pars$c_nugget,
                               pars$rho, n.c, n.c+n.e, 2) 




#### compare the fits ####

## score & MSE
#mse (probit)
MSE.hetgp <- MSE(y.test, Mean.het)
MSE.sml <- MSE(y.test, Mean.ml)
MSE.hetgp; MSE.sml

#mse (original scale)
MSE.hetgp0 <- MSE(y.test1, inv.probit(Mean.het))
MSE.sml0 <- MSE(y.test1, inv.probit(Mean.ml))
MSE.hetgp0; MSE.sml0

#score
Score.hetgp <- Score(y.test, Mean.het, diag(Var.het))
Score.sml <- Score(y.test, Mean.ml, diag(Var.ml))
sum(Score.hetgp); sum(Score.sml)
## look at residuals
par(mfrow=c(1,2))
resids.het <- (y.test - Mean.het)/sqrt(diag(Var.het))
resids.ml <- (y.test - Mean.ml)/sqrt(diag(Var.ml))

chol.het <- forwardsolve(t(chol(Var.het)), (y.test - Mean.het))
chol.ml <- forwardsolve(t(chol(Var.ml)), (y.test - Mean.ml))
hist(resids.het,prob=T)
hist(resids.ml,prob=T)

## coverage plot


coverage <- function(resids){
  alpha <- seq(0, 1, length = length(resids))
  emp.cov <- rep(0, length(resids))
  for(i in 1:length(resids)){
    emp.cov[i] <- sum(abs(resids) < qnorm(1-alpha[i]/2))/length(resids)
  }
  list(x=emp.cov,y = 1 - alpha) 
}
plot(y.test, Mean.ml)
par(mfrow = c(1,2))

plot(coverage(resids.het), xlab = "Empirical Coverage", ylab="", main = "HetGP", pch = 20)
abline(0,1)
title(ylab = "Expected Covereage" , line = 1.9)
plot(coverage(resids.ml),  xlab = "Empirical Coverage", ylab="", main = "SML", pch = 20)
abline(0,1)
par(mfrow=c(1,2), cex.axis=2, cex.lab=2,cex.main=2)
plot(coverage(chol.het), xlab = "Empirical Coverage", ylab="", main = "HetGP", pch = 20)
abline(0,1)
title(ylab = "Expected Covereage" , line = 2.4)
plot(coverage(chol.ml),  xlab = "Empirical Coverage", ylab="", main = "SML", pch = 20)
abline(0,1)