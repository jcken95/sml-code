## load stuff

## also fit the emulators!

#### read code ####
source("GPfunctionsOptim.R")
source("hetGPfunctions.R")
library(truncnorm)
probit <- function(p){
  qnorm(p)
}
inv.probit <- function(q){
  pnorm(q)
}

#### read data ####
pars <- readRDS("sml-params.RDS")
x.c <- as.matrix(read.table(file = "emulator-data/x-c.txt",sep=","))
x.e <- as.matrix(read.table(file = "emulator-data/x-e.txt",sep=","))
y.c1 <- read.table(file = "emulator-data/y-c.txt",sep=",")[,1]
y.e1 <- read.table(file = "emulator-data/y-e.txt",sep=",")[,1]

y.c <- probit(y.c1)
y.e <- probit(y.e1)


n.e <- dim(x.e)[1]
n.c <- dim(x.c)[1]
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))


set.seed(123) ## reproduceable analysis
n <- 15^2 ## number to sample in the x-direction
nDims <- 9
randomSample <- function(nn = n, nnDims = nDims){
  ## simulate RVs from expert's prior
  #matrix(runif(nn*nnDims,0.1,5),ncol=nnDims)
 lhs::randomLHS(nn,nnDims)*(5 - 0.1) + 0.1
}
xp <- randomSample()  # used for plotting techniques
saveRDS(xp, "xp.RDS")
N <- 15^2 ## number of Y_m(x) to esimate via vectorisation
## can be RAM heavy so don't make too big
## stick to max 10^4


mu.yobs <- c(cbind(1,log(x.c[n.c:1,]))%*%pars$c_beta, cbind(1,log(x.e[n.e:1,]))%*%(pars$rho*pars$c_beta + pars$m_beta))

mu.yobs.var <- cbind(1, x.e.std)%*%pars$v_beta
print("estimating f_0, Var(Y)")
y0 <- 0
mu.var <- as.vector(cbind(1,x.e.std)%*%pars$v_beta)
n.reps <- 100
for(i in 1:n.reps){
  
  x <- randomSample()
  x.std <- scale(x, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
  Var.sml <- cond.var.het(pars$logLambda, x.e.std, x.std,
                          pars$v_sigma, pars$v_theta, pars$v_nugget^2,
                          as.vector(cbind(1,x.std)%*%pars$v_beta), 
                          mu.var, 2)
  #Var.het <- diag( cond.covmat.het(exp(pars$logLambda), rep(0, 1000),  x.e.std, x,pars$v_sigma, pars$v_theta, 2) )
  y0 <- y0 + mean(Var.sml)
  print(i)
}
y0 <- y0/n.reps ## E(var(Y|X))
y0 ## this is expected variance

rm(Var.sml)
print("estimtated y_0")
## now find Var(E(Y|X))

y1 <- 0
f0 <- 0
n.reps <- 100

for(i in 1:n.reps){
  x <- randomSample()
  mu.y <-cbind(1,log(x))%*%(pars$rho*pars$c_beta + pars$m_beta)
  x.std <-  scale(x, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
  Mean.sml <- cond.mean.hetSML(c(y.c, y.e), rbind(x.c.std, x.e.std), x.std,
                   c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                   c(pars$m_sigma, pars$m_theta, 0), exp(pars$logLambda), pars$rho,
                   mu.y, mu.yobs, n.c, 2)
  
  y1 <- y1 + var(Mean.sml)
  f0 <- f0 + mean(Mean.sml)
  print(i)
}
y1 <- as.vector(y1/n.reps) ## Var(E(Y|x))
f0 <- as.vector(f0/n.reps) ## E(Y)
Var.y <- y0 + y1
100*y0/Var.y
print("estimtated Var(Y)")
#rm( Mean.het)
n.repeats <- 10
psa <- function(k,n,n.repeats,xp){
  
  means <- vars <-  rep(0,n)
  for(j in 1:n){ ## how many points to estimate variance?
    mean.temp <- var.temp <- 0
    
    for(l in 1:n.repeats){ ## how many points to estimate each expectation?
      
      x <- randomSample(nn = N)
      
      x[,k] <- rep(xp[j], N)
      
      mu.y <- cbind(1,log(x))%*%(pars$rho*pars$c_beta + pars$m_beta)
      x.std <- scale(x, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
      
      
      #mu.y.var <- cbind(1, x)%*%pars$v_beta
      
      
      y.temp <- rep(0, N)
      
      y.temp <- cond.mean.hetSML(c(y.c, y.e), rbind(x.c.std, x.e.std), x.std,
                                 c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                                 c(pars$m_sigma, pars$m_theta, 0), exp(pars$logLambda), pars$rho,
                                 mu.y, mu.yobs, n.c, 2)
      #var.temp0 <- mean(diag(cond.covmat.het(exp(pars$logLambda), rep(0,N),
      #                                       x.e.std, x, pars$m_sigma, pars$m_theta, nu=2)))
      mean.temp <-  mean(y.temp) + mean.temp
      var.temp <- var(y.temp)*(N-1)/N + var.temp
    }
    
    means[j] <- mean.temp/n.repeats
    vars[j] <- var.temp/n.repeats
    
  }
  res <- list(means = means, vars = vars,x = xp)
  res
}

psa.wrapper <- function(k, ...){
  psa(k, n, n.repeats, xp[,k])
}
st <- Sys.time()
psa.list <- parallel::mclapply(X = 1:9, FUN = psa.wrapper, n = 100, n.repeats = 1, xp = xp,mc.cores = 8)
en <- Sys.time()
en-st
psa.list
par(mfrow=c(1,1))
for(i in 1:9){
  plot(psa.list[[i]]$x, psa.list[[i]]$means, ylim = probit(c(0.75,0.98)),type="l")
}



## total uncertainty

saveRDS(psa.list,"psa-sml.RDS")
extras <- list(E.var.y = y0, Var.E.y = y1, E.y = f0)
saveRDS(extras, "extras.RDS")

