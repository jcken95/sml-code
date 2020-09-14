## compute the second order terms; f_{ij} (x_{ij}) = f(X|X_{ij}) - f_i(x_i) - f_j(x_j) - f_0
## just compute f(X|X_{ij}) in this term

## load stuff

## also fit the emulators!

#### read code ####
source("GPfunctionsOptim.R")
source("hetGPfunctions.R")
library(truncnorm)
library("compiler")
#### read data ####
pars <- readRDS("sml-params.RDS")
x.c <- as.matrix(read.table(file = "emulator-data/x-c.txt",sep=","))
x.e <- as.matrix(read.table(file = "emulator-data/x-e.txt",sep=","))
y.c1 <- read.table(file = "emulator-data/y-c.txt",sep=",")[,1]
y.e1 <- read.table(file = "emulator-data/y-e.txt",sep=",")[,1]

probit <- function(p) qnorm(p)

y.c <- probit(y.c1)
y.e <- probit(y.e1)


n.e <- dim(x.e)[1]
n.c <- dim(x.c)[1]
x.e.std <- scale(x.e)
x.c.std <- scale(x.c, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
xp <- readRDS("xp.RDS")
#comp.mean <- cmpfun(cond.mean.hetSML)
set.seed(1234) ## reproduceable analysis
n <- 15^2 ## number to sample in the x-direction
nDims <- 9
randomSample <- function(nn = n, nnDims = nDims){
  ## simulate RVs from expert's prior
  #matrix(runif(nn*nnDims,0.1,5),ncol=nnDims)
  matrix(rtruncnorm(nn*nnDims, 0.1,5,2.51, 1.205),ncol=nnDims)
}
#saveRDS(xp, "xp.RDS")
N <- 15^2 ## number of Y_m(x) to esimate via vectorisation
## can be RAM heavy so don't make too big
## stick to max 10^4

#xp <- randomSample()
mu.yobs <- c(cbind(1,log(x.c[n.c:1,]))%*%pars$c_beta, cbind(1,log(x.e[n.e:1,]))%*%(pars$rho*pars$c_beta + pars$m_beta))

mu.yobs.var <- cbind(1, x.e.std)%*%pars$v_beta
print("estimating f_0, Var(Y)")
y0 <- 0
mu.var <- as.vector(cbind(1,x.e.std)%*%pars$v_beta)

psa <- function(i,N,n.repeats,xp, K=9){
  means <- list()
  vars <- list()
  xi <- xp[,i]
  for(j in (i+1):K){
    means.temp <- vars.temp <- matrix(0, nrow = N, ncol = N)
    xj <- xp[,j]
    for(reps in 1:n.repeats){
      mean.temp <- var.temp <- 0
      for(ii in 1:N){
        for(jj in 1:N){

          x <- randomSample()
          x[,i] <- rep(xi[ii],N)
          x[,j] <- rep(xj[jj],N)
          mu.y <- cbind(1,log(x))%*%(pars$rho*pars$c_beta + pars$m_beta)
          x.std <- scale(x, center = attr(x.e.std, "scaled:center"), scale = attr(x.e.std, "scaled:scale"))
          
          y.temp <- rep(0, N)
          y.temp <- cond.mean.hetSML(c(y.c, y.e), rbind(x.c.std, x.e.std), x.std,
                                     c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                                     c(pars$m_sigma, pars$m_theta, 0), exp(pars$logLambda), pars$rho,
                                     mu.y, mu.yobs, n.c, 2)

          mean.temp <-  mean(y.temp)
          var.temp <- var(y.temp)
          means.temp[ii,jj] <- means.temp[ii,jj] + mean.temp/n.repeats
          vars.temp[ii,jj] <- vars.temp[ii,jj] + var.temp/n.repeats
          
        }
        

      }

    }
  means[[j]] <- means.temp
  vars[[j]] <- vars.temp
  }
  

  
  res <- list(means = means, vars = vars,x = xp)
  res
}

psa.wrapper <- function(i, ...){
  psa(i,n,n.repeats,xp, K=9)
}
st <- Sys.time()
interaction.list <- parallel::mclapply(X = 1:8, FUN = psa.wrapper, n = 10, n.repeats = 1, xp = xp,mc.cores = 8)
en <- Sys.time()
en-st
interaction.list
par(mfrow=c(1,1))
for(i in 1:9){
  plot(psa.list[[i]]$x, psa.list[[i]]$means, ylim = probit(c(0.75,0.98)),type="l")
}



## total uncertainty

saveRDS(interaction.list,"interactions.RDS")
#extras <- list(E.var.y = y0, Var.E.y = y1, E.y = f0)
#saveRDS(extras, "sml-sa/extras.RDS")

