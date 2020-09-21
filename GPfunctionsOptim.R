## R functions for Gaussian Process based emulation
## 11/6/10


## MVN functions
#function(psi,psi.lower,psi.upper)
#{
  ## uniform priors
  
#  -sum(log(psi.upper[-1]-psi.lower[-1]))
#}

matfun <- function(A,fun)
{
  ## applies a function fun (which expects a single scalar argument)
  ## to a matrix A, in a sensible way...
  s <- try(svd(A))
  ##s <- svd(A)
  ( s$u %*% diag( fun(s$d) ) %*% s$v )
}

mat.sqrt.dh <- function(A)
{
  
  ## computes the matrix square root of the symmetric positive
  ## (semi-)definite matrix A
  
  s <- try(svd(A))
  
  if(!is.numeric(s$d)){
    
    print("spectral decomposition failed")
    sqrt.A <- NULL
  }
  else{
    sqrt.A <-  s$u %*% diag( sqrt(s$d) ) %*% s$v 
  }
  
  if(is.null(sqrt.A)) sqrt.A <- chol(A)
  
  sqrt.A
}

mat.sqrt <- function(A){
  chol(A)
}



logdmvnorm.sqrt2 <- function(x,mu,sqrt.Sigma) ## my way
{
  ## Compute the log density of MVN random vector x with mean mu
  ## and covariance matrix sqrt.Sigma%*%sqrt.Sigma.
  ## Faster than corresponding function logdmvnorm.
  
  ## NB expects that A is upper/lower triang. (eg cholesky)
  
  ##
  ## created: 30/08/06
  ## updated: 28/11/06
  ## updated: 20/12/18 (JK)
  
  n <- length(x)
  A <- sqrt.Sigma
  
  log.det.A <- sum(log(diag(A)))
  
  b <- solve(A,(x-mu))
  bb <- t(b)%*%b
  as.numeric(-0.5*(n*log(2*pi)+2*log.det.A+bb))
}


logdmvnorm.sqrt <- function(x,mu,sqrt.Sigma) ## my way
{
  ## Compute the log density of MVN random vector x with mean mu
  ## and covariance matrix sqrt.Sigma%*%sqrt.Sigma.
  ## Faster than corresponding function logdmvnorm.
  
  ## NB expects that A is upper/lower triang. (eg cholesky)
  
  ##
  ## created: 30/08/06
  ## updated: 28/11/06
  ## updated: 20/12/18 (JK)
  ## updated: 4/1/19 (JK)
  
  n <- length(x)
  
  
  log.det.A <- sum(log(diag(sqrt.sigma)))
  
  b <- backsolve(A,(x-mu))
  bb <- t(b)%*%b
  as.numeric(-0.5*(n*log(2*pi)+2*log.det.A+bb))
}



rmvnorm.sqrt <- function(n,mu,sqrt.Sigma)
{
  ## Generate a random sample of size n from the multivariate normal
  ## distribution with mean vector mu and covariance matrix
  ## sqrt.Sigma%*%sqrt.Sigma
  ## Faster than corresponding rmvnorm in GPbasics.
  ##
  ## created: 30/08/06
  ## updated: 28/11/06
  
  p <- dim(sqrt.Sigma)[1]
  eps <- matrix(rnorm(n*p),nrow=p,ncol=n)
  mu %*% t(rep(1,n)) + ( sqrt.Sigma %*% eps )
}


cor.kern.fun <- 
  function(x1,x2,w=1,p=2)
  {
    ## The 'kernel' of the correlation function
    ## w is a positive roughness parameter
    ## p is a positve power (0 < p <= 2)  
    
    (abs((x1-x2)/w)^p)
  }



cov.x1.x2 <- 
  function(x1,x2,sigma=1,omega=rep(1,dim(x1)[2]),nu=2)
  { 
    ## generates matrix of covariances between coordinates x1 and x2
    ## for example, x1 might be the matrix containing coordinates of
    ## the points at which prediction is required
    ## and x2 might be the coordinates of the points at which data is observed
    
    dims.x1 <- dim(as.matrix(x1))
    mat <- matrix(0,nrow=dims.x1[1],ncol=dim(as.matrix(x2))[1])
    for(i in 1:dims.x1[2]){
      mat <- mat+outer(as.matrix(x1)[,i],as.matrix(x2)[,i],cor.kern.fun,w=omega[i],p=nu)
    }
    (sigma^2)*exp(-mat)
  }

log.cov.x1.x2 <- function(x1,x2,sigma=1,omega=rep(1,dim(x1)[2]),nu=2){
  dims.x1 <- dim(as.matrix(x1))
  mat <- matrix(0,nrow=dims.x1[1],ncol=dim(as.matrix(x2))[1])
  for(i in 1:dims.x1[2]){
    mat <- mat+outer(as.matrix(x1)[,i],as.matrix(x2)[,i],cor.kern.fun,w=omega[i],p=nu)
  }
  
  -mat
  
}


cond.mean <- function(y.obs,x.obs,x.new,sigma,omega,mu.Y,mu.Yobs,nugget=0,nu=2)
{
  
  ## computes the conditional mean at points x.new
  ## based on observations y.obs at points x.obs  
  
  ## assuming independent measurement error with variance = nugget^2 
  
  Sigma.Yobs.Yobs <- cov.x1.x2(x.obs,x.obs,sigma,omega,nu)
  Sigma.Y.Yobs <- cov.x1.x2(x.new,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs+diag(nugget^2,dim(Sigma.Yobs.Yobs)[1])
  
  ## find solve(Sigma.Yobs.Yobs) carefully...
  
  R <- chol(Sigma.Yobs.Yobs) ##sqrt s.yobs.yobs
  
  R.inv <- backsolve(R, diag(1, dim(R)[1]))
  
  mu.Y+Sigma.Y.Yobs %*%  ((R.inv) %*% (t(R.inv) %*%(y.obs-mu.Yobs)))
}  

cond.covmat <- function(y.obs,x.obs,x.new,sigma,omega,nugget=1e-10,nu=2)
{
  ## computes the conditional covariance matrix at points x.new
  ## based on observations y.obs at points x.obs  
  
  ## assuming independent measurement error with variance = nugget 
  
  nug.sq <- nugget^2
  Sigma.Y.Y <- cov.x1.x2(x.new,x.new,sigma,omega,nu)
  Sigma.Y.Y <- Sigma.Y.Y+diag(nug.sq,dim(Sigma.Y.Y)[1])
  Sigma.Y.Yobs <- cov.x1.x2(x.new,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- cov.x1.x2(x.obs,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs+diag(nug.sq,dim(Sigma.Yobs.Yobs)[1])
  
  R <- chol(Sigma.Yobs.Yobs) ##sqrt s.yobs.yobs
  
  R.inv <- backsolve(R, diag(1, dim(R)[1])) ## inv sqrt s.yobs.yobs
  
  
  Sigma.Y.Y- ((Sigma.Y.Yobs) %*% (R.inv))  %*% t((Sigma.Y.Yobs) %*% (R.inv)) ## syy - s.y.yobs * inv s.yobs.yobs * s.yobs.y
}  


cond.mean.old <- function(y.obs,x.obs,x.new,sigma,omega,mu.Y,mu.Yobs,nugget=0,nu=2)
{
  ## old version - don't use
  ## computes the conditional mean at points x.new
  ## based on observations y.obs at points x.obs  
  
  ## assuming independent measurement error with variance = nugget^2 
  
  Sigma.Yobs.Yobs <- cov.x1.x2(x.obs,x.obs,sigma,omega,nu)
  Sigma.Y.Yobs <- cov.x1.x2(x.new,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs+diag(nugget^2,dim(Sigma.Yobs.Yobs)[1])
  if(dim(Sigma.Yobs.Yobs)[1]==1){
    sq <- sqrt(Sigma.Yobs.Yobs[1,1])
  }
  else{
    sq <- mat.sqrt(Sigma.Yobs.Yobs)
    #matfun(Sigma.Yobs.Yobs,sqrt)
  }
  a <- t(solve(sq,t(Sigma.Y.Yobs)))
  ##  a <- Sigma.Y.Yobs%*%solve(t(sq))
  b <- backsolve(sq,(y.obs-mu.Yobs))
  mu.Y+a%*%b
}  


cond.covmat.old <- function(y.obs,x.obs,x.new,sigma,omega,nugget=1e-10,nu=2)
{
  ## old version - don't use - introduces errors!
  ## computes the conditional covariance matrix at points x.new
  ## based on observations y.obs at points x.obs  
  
  ## assuming independent measurement error with variance = nugget 
  
  nug.sq <- nugget^2
  Sigma.Y.Y <- cov.x1.x2(x.new,x.new,sigma,omega,nu)
  Sigma.Y.Y <- Sigma.Y.Y++diag(nug.sq,dim(Sigma.Y.Y)[1])
  Sigma.Yobs.Yobs <- cov.x1.x2(x.obs,x.obs,sigma,omega,nu)
  Sigma.Y.Yobs <- cov.x1.x2(x.new,x.obs,sigma,omega,nu)
  Sigma.Yobs.Yobs <- Sigma.Yobs.Yobs+diag(nug.sq,dim(Sigma.Yobs.Yobs)[1])
  if(dim(Sigma.Yobs.Yobs)[1]==1){
    sq <- sqrt(Sigma.Yobs.Yobs[1,1])
  }
  else{
    sq <- matfun(Sigma.Yobs.Yobs,sqrt)
  }
  b <- solve(sq,t(Sigma.Y.Yobs))
  Sigma.Y.Y-t(b)%*%b
}  
?solve

neg.ll <-  function(psi,y,x,nu=2)
{
  n.psi <- length(psi)
  mu <- psi[1]
  sigma <- exp(psi[2])
  omega <- exp(psi[-c(1,2,n.psi)])
  lambda <- exp(psi[n.psi])
  
  nll <- -llike.GP(y,x,mu,sigma,omega,lambda,nu) 
  print(c(psi,nll))
  nll
  
}

neg.ll1 <-  function(psi,y,x,nu=2,mu=0)
{
  n.psi <- length(psi)
  sigma <- exp(psi[1])
  omega <- exp(psi[-c(1,n.psi)])
  lambda <- exp(psi[n.psi])
  
  nll <- -llike.GP(y,x,mu,sigma,omega,lambda,nu) 
  print(c(psi,nll))
  nll
  
}

neg.ll2 <-  function(psi,y,x,nu=2)
{
  n.psi <- length(psi)
  omega <- exp(psi[-n.psi])
  lambda <- exp(psi[n.psi])
  
  nll <- -llike.GP(y,x,mu=0,sigma=sd(y),omega,lambda,nu) 
  print(c(psi,nll))
  nll
  
}


inv.part <- function(A,B,C.inv)
{ 
  ## inverse of a partitioned matrix
  B.t <- t(B)
  C.inv.mult.B <- C.inv%*%B
  Q <- A-B.t%*%(C.inv.mult.B)
  Q.inv <- solve(as.matrix(Q))
  A.up <- Q.inv
  B.up <- -C.inv.mult.B%*%Q.inv
  B.up.t <- t(B.up)
  C.up <- C.inv+(-B.up%*%t(C.inv.mult.B))
  X.inv <- cbind(rbind(A.up,B.up),rbind(B.up.t,C.up))
  X.inv
}

inv.part.no.inv <- function(A,B,C)
{ 
  ## inverse of a partitioned matrix
  C.inv <- solve(C)
  B.t <- t(B)
  C.inv.mult.B <- C.inv%*%B
  Q <- A-B.t%*%(C.inv.mult.B)
  Q.inv <- solve(as.matrix(Q))
  A.up <- Q.inv
  B.up <- -C.inv.mult.B%*%Q.inv
  B.up.t <- t(B.up)
  C.up <- C.inv+(-B.up%*%t(C.inv.mult.B))
  X.inv <- cbind(rbind(A.up,B.up),rbind(B.up.t,C.up))
  X.inv
}

dlexp.log <- function(y,theta)
{
  ## log of density function of Y=log X when X~Exp(theta)
  
  log(theta)+(y-theta*exp(y))
}


lprior.exp <- function(psi,a.psi)
{
  ## exponential priors
  
  sum(dlexp.log(psi[-1],a.psi[-1]))
}

lprior.unif.psi <- function(psi,psi.lower,psi.upper)
{
  ## uniform priors
  
  -sum(log(psi.upper[-1]-psi.lower[-1]))
}

lprior.unif <- function(theta,theta.lower,theta.upper)
{
  ## uniform priors
  
  -sum(log(theta.upper-theta.lower))
}

lprior.norm <- function(b.mean,b.cov,beta){
  sig.inv = solve(b.cov)
  -0.5*determinant(2*pi*b.cov) - 0.5*t(beta-b.mean)%*%sig.inv%*%(beta-b.mean)
}

llike.GP = function(y,x,mu,sigma,omega,lambda,nu=2)
{
  
  ## Compute the loglikelihood of parameters in a GP based on data y
  ## at x.
  
  ## created: 30/08/06
  ## updated: 1/10/09
  
  n <- length(y)
  Sigma.Yfull.Yfull <- cov.x1.x2(x,x,sigma,omega,nu)+diag(lambda^2,n)
  if(n==1){
    A <- sqrt(Sigma.Yfull.Yfull)
  }
  else {
    ##    A <- matfun(Sigma.Yfull.Yfull,sqrt)
    A <- mat.sqrt(Sigma.Yfull.Yfull) ## makes det calc easier
  }
  if(is.numeric(A)){
    if(n == 1){  
      log.det.A <- determinant(A,logarithm=TRUE)$modulus[1]
      b <- solve(A,(y-mu))
      bb <- t(b)%*%b
      ##  print(c(-(0.5*n*log(2*pi)),-log.det.X,-0.5*bb))
      res <- as.vector(-0.5*(n*log(2*pi)+2*log.det.A+bb))
    }else{
      
      log.det.A <- abs(sum(log(diag(A))))
      b <- backsolve(A,(y-mu))
      bb <- t(b)%*%b
      ##  print(c(-(0.5*n*log(2*pi)),-log.det.X,-0.5*bb))
      res <- as.vector(-0.5*(n*log(2*pi)+2*log.det.A+bb))
    }
    
  }
  else{
    ## the spectral decomposition failed and returned NULL
    res <- NULL
  }
  res
}


#llike.GP = function(y,x,mu = 0,sigma,omega,lambda,nu=2, trend = F, beta0 = 0, beta1 = 0)
#  {

## Compute the loglikelihood of parameters in a GP based on data y
## at x.

## created: 30/08/06
## updated: 2/10/18 (JK)
## allows for linear trend function
#    
#    n <- length(y)
#    if(trend){ mu <- beta0 + beta1*x }
#   
#    Sigma.Yfull.Yfull <- cov.x1.x2(x,x,sigma,omega,nu)+diag(lambda^2,n)
#    if(n==1){
#      A <- sqrt(Sigma.Yfull.Yfull)
#    }
#    else {
##    A <- matfun(Sigma.Yfull.Yfull,sqrt)
#      A <- mat.sqrt(Sigma.Yfull.Yfull)
#    }
#    if(is.numeric(A)){
#      log.det.A <- determinant(A,logarithm=TRUE)$modulus[1]
#      b <- solve(A,(y-mu))
#      bb <- t(b)%*%b
##  print(c(-(0.5*n*log(2*pi)),-log.det.X,-0.5*bb))
#      res <- as.vector(-0.5*(n*log(2*pi)+2*log.det.A+bb))
#    }
#    else{
## the spectral decomposition failed and returned NULL
#      res <- NULL
#    }
#    res
#  }

mcmc.GP <- function(y,x,nu=2,psi.curr,psi.lower,psi.upper,sqrt.prop.cov.mat,p.prior,its=100,burn=0,thin=1)
{
  ##I THINK that nu = 2 gives squared exp. kernel
  ## MCMC algorithm for sampling from posterior
  ## distribution of log-transformed GP hyperparameters
  ## assuming a unifrom prior
  
  ## the number of hyperparameters (mu,sigma,omega,lambda)
  n.psi <- length(psi.curr)
  
  ## current log-likelihood and log prior density
  llike.curr <- llike.GP(y,x,psi.curr[1],exp(psi.curr[2]),exp(psi.curr[-c(1,2,n.psi)]),exp(psi.curr[n.psi]),nu)
  lprior.curr  <- lprior.unif.psi(psi.curr,psi.lower,psi.upper)
  
  
  ## set up matrices/vectors to store the results
  psi.res <- matrix(0,nrow=its,ncol=n.psi)
  llike.res <- rep(0,its)
  lprior.res <- rep(0,its)
  
  for(i in 1:(burn+its)){
    
    print(i-burn)
    
    for(j in 1:thin){
      
      if(runif(1,0,1) <= p.prior){
        ## proposal from the prior
        
        psi.prop <- runif(n.psi,psi.lower,psi.upper)
        
        llike.prop <- llike.GP(y,x,psi.prop[1],exp(psi.prop[2]),exp(psi.prop[-c(1,2,n.psi)]),exp(psi.prop[n.psi]),nu)
        
        lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper)
        
        l.acc <- llike.prop-llike.curr
        
        if(log(runif(1,0,1)) < min(l.acc,0)){
          
          psi.curr <- psi.prop
          llike.curr <- llike.prop
          lprior.curr <- lprior.prop
        }
      }
      else{
        ## random walk proposal 
        
        psi.prop <- rmvnorm.sqrt(1,psi.curr,sqrt.prop.cov.mat)
        
        ##prop.is.valid <- prod((psi.prop < psi.upper)[-1]*(psi.prop > psi.lower)[-1])
        prop.is.valid <- prod((psi.prop < psi.upper)*(psi.prop > psi.lower)) 
        #product of logicals indicating whether proposal inside the parameter space
        if(prop.is.valid){
          
          llike.prop <- llike.GP(y,x,psi.prop[1],exp(psi.prop[2]),exp(psi.prop[-c(1,2,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper)
          
          l.acc <- lprior.prop-lprior.curr+llike.prop-llike.curr
          
          if(log(runif(1,0,1)) < min(l.acc,0)){
            
            psi.curr <- psi.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
          }
        }
      }
    }
    
    if(i > burn){
      ind <- i-burn
      psi.res[ind,] <- psi.curr
      llike.res[ind] <- llike.curr
      lprior.res[ind] <- lprior.curr
    }
  }
  list(psi=psi.res,llike=llike.res,lprior=lprior.res)
}

powerpost.GP <- function(y,x,nu=2,psi.curr,psi.lower,psi.upper,sqrt.prop.cov.mat,p.prior,its=100,burn=0,thin=1,temp=1)
{
  
  ## MCMC algorithm for sampling from power posterior
  ## distribution of log-transformed GP hyperparameters
  ## assuming a unifrom prior
  
  ## the number of hyperparameters (mu,sigma,omega,lambda)
  n.psi <- length(psi.curr)
  
  ## current log-likelihood and log prior density
  llike.curr <- llike.GP(y,x,psi.curr[1],exp(psi.curr[2]),exp(psi.curr[-c(1,2,n.psi)]),exp(psi.curr[n.psi]),nu)
  lprior.curr  <- lprior.unif.psi(psi.curr,psi.lower,psi.upper)
  
  
  ## set up matrices/vectors to store the results
  psi.res <- matrix(0,nrow=its,ncol=n.psi)
  llike.res <- rep(0,its)
  lprior.res <- rep(0,its)
  
  for(i in 1:(burn+its)){
    
    print(i-burn)
    
    for(j in 1:thin){
      
      if(runif(1,0,1) <= p.prior){
        ## proposal from the prior
        
        psi.prop <- runif(n.psi,psi.lower,psi.upper)
        
        llike.prop <- llike.GP(y,x,psi.prop[1],exp(psi.prop[2]),exp(psi.prop[-c(1,2,n.psi)]),exp(psi.prop[n.psi]),nu)
        
        lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper)
        
        l.acc <- temp*(llike.prop-llike.curr)
        
        if(log(runif(1,0,1)) < min(l.acc,0)){
          
          psi.curr <- psi.prop
          llike.curr <- llike.prop
          lprior.curr <- lprior.prop
        }
      }
      else{
        ## random walk proposal 
        
        psi.prop <- rmvnorm.sqrt(1,psi.curr,sqrt.prop.cov.mat)
        
        prop.is.valid <- prod((psi.prop < psi.upper)[-1]*(psi.prop > psi.lower)[-1])
        if(prop.is.valid){
          
          llike.prop <- llike.GP(y,x,psi.prop[1],exp(psi.prop[2]),exp(psi.prop[-c(1,2,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper)
          
          l.acc <- lprior.prop-lprior.curr+temp*(llike.prop-llike.curr)
          
          if(log(runif(1,0,1)) < min(l.acc,0)){
            
            psi.curr <- psi.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
          }
        }
      }
    }
    
    if(i > burn){
      ind <- i-burn
      psi.res[ind,] <- psi.curr
      llike.res[ind] <- llike.curr
      lprior.res[ind] <- lprior.curr
    }
  }
  list(psi=psi.res,llike=llike.res,lprior=lprior.res)
}

parallel.temp.GP <- function(y,x,nu=2,psi.curr,psi.lower,psi.upper,sqrt.prop.cov.mat,p.prior,its=100,burn=0,thin=1,temp=1,swap=TRUE)
{
  
  ## parallel tempering algorithm for sampling from posterior
  ## distribution of log-transformed GP hyperparameters
  ## assuming a unifrom prior
  
  ## number of chains
  M <- length(temp)
  
  ## indices of chains to swap
  ind.swap <- rep(0,2)
  
  ## the number of hyperparameters (mu,sigma,omega,lambda)
  n.psi <- dim(psi.curr)[2]
  
  llike.curr <- rep(0,M)
  lprior.curr <- rep(0,M)  
  for(m in 1:M){
    ## current log-likelihood and log prior density
    llike.curr[m] <- llike.GP(y,x,psi.curr[m,1],exp(psi.curr[m,2]),exp(psi.curr[m,-c(1,2,n.psi)]),exp(psi.curr[m,n.psi]),nu)
    lprior.curr[m]  <- lprior.unif.psi(psi.curr[m,],psi.lower,psi.upper)
  }
  
  ## set up matrices/vectors to store the results
  psi.res <- array(0,c(its,n.psi,m))
  llike.res <- matrix(0,nrow=its,ncol=m)
  lprior.res <- matrix(0,nrow=its,ncol=m)
  
  for(i in 1:(burn+its)){
    
    print(i-burn)
    
    for(j in 1:thin){
      
      
      for(m in 1:M){
        
        if(runif(1,0,1) <= p.prior[m]){
          ## proposal from the prior
          
          psi.prop <- runif(n.psi,psi.lower,psi.upper)
          
          llike.prop <- llike.GP(y,x,psi.prop[1],exp(psi.prop[2]),exp(psi.prop[-c(1,2,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper)
          
          l.acc <- temp[m]*(llike.prop-llike.curr[m])
          
          if(log(runif(1,0,1)) < min(l.acc,0)){
            
            psi.curr[m,] <- psi.prop
            llike.curr[m] <- llike.prop
            lprior.curr[m] <- lprior.prop
          }
        }
        else{
          ## random walk proposal 
          
          psi.prop <- rmvnorm.sqrt(1,psi.curr[m,],sqrt.prop.cov.mat[,,m])
          
          prop.is.valid <- prod((psi.prop < psi.upper)[-1]*(psi.prop > psi.lower)[-1])
          if(prop.is.valid){
            
            llike.prop <- llike.GP(y,x,psi.prop[1],exp(psi.prop[2]),exp(psi.prop[-c(1,2,n.psi)]),exp(psi.prop[n.psi]),nu)
            
            lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper)
            
            l.acc <- lprior.prop-lprior.curr[m]+temp[m]*(llike.prop-llike.curr[m])
            
            if(log(runif(1,0,1)) < min(l.acc,0)){
              
              psi.curr[m,] <- psi.prop
              llike.curr[m] <- llike.prop
              lprior.curr[m] <- lprior.prop
            }
          }
        }
      }
    }
    
    if(swap){
      if(M>1){
        ind.swap[1] <- sample(1:(M-1),1,replace=FALSE)
        ind.swap[2] <- ind.swap[1]+1
        
        l.acc.swap <- llike.curr[ind.swap[1]]*(temp[ind.swap[2]]-temp[ind.swap[1]])+llike.curr[ind.swap[2]]*(temp[ind.swap[1]]-temp[ind.swap[2]])
        
        if(log(runif(1,0,1)) < min(l.acc.swap,0)){
          print(paste("swap",ind.swap))
          
          psi.new <- psi.curr[ind.swap[1],]
          psi.curr[ind.swap[1],] <- psi.curr[ind.swap[2],]
          psi.curr[ind.swap[2],] <- psi.new
          
          llike.new <- llike.curr[ind.swap[1]]
          llike.curr[ind.swap[1]] <- llike.curr[ind.swap[2]]
          llike.curr[ind.swap[2]] <- llike.new
          
          lprior.new <- lprior.curr[ind.swap[1]]
          lprior.curr[ind.swap[1]] <- lprior.curr[ind.swap[2]]
          lprior.curr[ind.swap[2]] <- lprior.new
          
        }
      }
    }
    
    if(i > burn){
      ind <- i-burn
      psi.res[ind,,] <- t(psi.curr)
      llike.res[ind,] <- llike.curr
      lprior.res[ind,] <- lprior.curr
    }
  }
  list(psi=psi.res,llike=llike.res,lprior=lprior.res)
}

###########
#JK's GP fn which includes a mean fn

mcmc.mean.gp <- 
  function(y, x, nu=2, psi.curr, psi.lower, psi.upper, beta.curr, cov.beta.prop.sqrt, sqrt.prop.cov.mat, beta.prior.mean, cov.beta.sqrt, p.prior, its=100, burn=0, thin=1){
    
    #psi = (log sig, log omega, log lambda) - log hyper param
    #beta = (b0,b1,...,bn) coeffs of basis fns
    ## mcmc scheme for sampling from posterior
    ## when the mean of the gp a priori is (potentially) non zero
    ## assuming a unifrom prior, but the beta are N(b_0, sig_beta)
    
    
    beta.prior.cov.sqrt <- cov.beta.sqrt
    n = length(y)
    ones = rep(1,n)
    X = cbind(ones,x) ##for now, we must have x as the transformed data (so, x := (h_0(x). h_1(x),...))
    n.psi <- length(psi.curr)
    n.beta <- length(beta.curr)
    mu.curr <- X %*% beta.curr ##might change this exact formatting
    
    ## current log-likelihood and log prior density
    
    llike.curr <- llike.GP(y,x,mu.curr,exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),nu)
    lprior.curr  <- lprior.unif.psi(psi.curr,psi.lower,psi.upper) + logdmvnorm.sqrt(beta.curr, mu = beta.curr, sqrt.Sigma = cov.beta.sqrt)
    
    
    ##store results in the following.sqrt
    psi.res <- matrix(0, nrow = its, ncol = n.psi)
    beta.res <- matrix(0, nrow = its, ncol = n.beta)
    llike.res <- rep(0, its)
    lprior.res <- rep(0, its)
    
    for(i in 1:(burn+its)){
      print(i-burn) #its counter
      
      for(j in 1:thin){
        if(runif(1,0,1) <= p.prior){
          ##propose from the prior
          psi.prop <- c(runif(n.psi,psi.lower,psi.upper))
          
          beta.prop <- rmvnorm.sqrt(1, mu = beta.prior.mean, sqrt.Sigma = cov.beta.sqrt)
          
          mu.prop <- X %*% beta.prop
          
          llike.prop <- llike.GP(y,x,mu.prop,exp(psi.prop[1]),exp(psi.prop[-c(1,n.psi)]),exp(psi.prop[n.psi]),nu)
          ##assume psi indep. beta a prior => log prior density is sum of log of components
          lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper) + logdmvnorm.sqrt(beta.prop, mu = beta.prior.mean, sqrt.Sigma = cov.beta.sqrt)
          
          
          l.acc <- llike.prop-llike.curr
          
          if(log(runif(1,0,1)) < min(l.acc,0)){
            ##will we move this way?
            psi.curr <- psi.prop
            beta.curr <- beta.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
            mu.curr <- X %*% beta.curr
          }
          
        }else{
          ## random walk proposal 
          
          psi.prop <- rmvnorm.sqrt(1, mu = psi.curr, sqrt.Sigma = sqrt.prop.cov.mat) #might exceed the param space, hence validate
          beta.prop <- rmvnorm.sqrt(1, mu =  beta.curr,  sqrt.Sigma = cov.beta.prop.sqrt) #always a valid prop
          mu.prop <- X %*% beta.prop #value of mu at the proposed unknowns
          
          ##prop.is.valid <- prod((psi.prop < psi.upper)[-1]*(psi.prop > psi.lower)[-1])
          prop.is.valid <- prod((psi.prop < psi.upper)*(psi.prop > psi.lower)) 
          #product of logicals indicating whether proposal inside the parameter space
          if(prop.is.valid){
            
            llike.prop <- llike.GP(y,x,mu.prop,exp(psi.prop[1]),exp(psi.prop[-c(1,n.psi)]),exp(psi.prop[n.psi]),nu)
            lprior.prop  <- lprior.unif.psi(psi.prop,psi.lower,psi.upper) + logdmvnorm.sqrt(beta.prop, mu = beta.curr, sqrt.Sigma = cov.beta.sqrt)
            
            
            l.acc <- lprior.prop-lprior.curr+llike.prop-llike.curr
            
            if(log(runif(1,0,1)) < min(l.acc,0)){
              
              psi.curr <- psi.prop
              beta.curr <- beta.prop
              llike.curr <- llike.prop
              lprior.curr <- lprior.prop
              
            }
          }
        }
      }
      if(i > burn){
        ind <- i-burn
        psi.res[ind,] <- psi.curr
        beta.res[ind,] <- beta.curr
        llike.res[ind] <-  llike.curr
        lprior.res[ind] <- lprior.curr
      }
    }
    list(psi = psi.res, beta = beta.res, llike=llike.res, lprior=lprior.res)
  }

## log normal mcmc scheme (hyper params are LN~(.,.))
## for an emulator with one input
## also this is a step-wise mcmc scheme (in psi)
mcmc.ln = function(y, x, nu=2, psi.curr, psi.mean, cov.psi.sqrt,
                   beta.curr, cov.beta.prop.sqrt, sqrt.prop.cov.mat,
                   beta.prior.mean, cov.beta.sqrt, p.prior, its=100, 
                   burn=0, thin=1){
  
  ## psi = (mu, log sig, log omega_vec, log lambda)
  beta.prior.cov.sqrt <- cov.beta.sqrt
  n <- length(y)
  ones <- rep(1,n)
  X <- cbind(ones,x) ##for now, we must have x as the transformed data (so, x := (h_0(x). h_1(x),...))
  
  n.psi <- length(psi.curr)
  n.beta <- length(beta.curr)
  mu.curr <- X %*% beta.curr #current value of the mean function at the design points
  
  
  llike.curr <- llike.GP(y,x,mu.curr,exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),nu)
  ## will assume log(hyper params) ~ N(.,.)
  lprior.curr <- logdmvnorm.sqrt(psi.curr,psi.mean,cov.psi.sqrt) + logdmvnorm.sqrt(beta.curr,beta.prior.mean,beta.prior.cov.sqrt)
  
  ## will work on log scale for hyper param
  
  logSig <- rep(0,its)
  n.Omega <- n.psi-2
  logOmega <- matrix(rep(0,its*n.Omega),ncol = n.Omega)
  logLambda <- rep(0,its)
  
  psi.res = cbind(logSig,logOmega,logLambda)
  
  ## beta will be on the standard scale
  
  n.beta = length(beta.curr)
  
  beta.res = matrix(rep(0,n.beta*its),ncol = n.beta, nrow = its)
  
  
  llike.res <- rep(0,its)
  lprior.res <- rep(0,its)
  
  for(i in 1:(burn+its)){
    
    print(i-burn)
    
    for(j in 1:thin){
      
      ## first work through each component of psi
      for(k in 1:n.psi){
        ## propose from prior
        
        if(runif(1,0,1) <= p.prior){
          psi.prop <- psi.curr
          
          psi.prop[k] <- rnorm(1,psi.mean[k],cov.psi.sqrt[k,k])
          
          llike.prop <- llike.GP(y,x,mu.curr,exp(psi.prop[1]),exp(psi.prop[-c(1,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          lprior.prop <- logdmvnorm.sqrt(psi.prop,psi.mean,cov.psi.sqrt)
          
          l.acc <- llike.prop - llike.curr 
          
          if(log(runif(1,0,1))  <= min(l.acc,0)){
            psi.curr <- psi.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
          }
          
        }else{
          ## random walk prop
          
          psi.prop <- psi.curr
          
          innovation <- rnorm(1,0,cov.psi.sqrt[k,k])
          
          psi.prop[k] <-  psi.curr[k] + innovation
          
          
          ##NB this is always a 'valid' proposal
          
          
          llike.prop <- llike.GP(y,x,mu.curr,exp(psi.prop[1]),exp(psi.prop[-c(1,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          ## now eval acc prob
          
          l.acc <- llike.prop - lprior.curr + llike.prop - llike.curr ## sym chain therefore do not consider q(.,.)
          
          if( log(runif(1,0,1)) <= min(l.acc,0)){
            psi.curr <- psi.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
          }
          
        }
        
      }
      ##now for the betas which will again be done in steps
      for(k in 1:n.beta){
        ## first propose from prior
        if(runif(1,0,1) <= p.prior){
          
          ## prop the k-th component of beta
          
          bbeta <- rnorm(1,beta.prior.mean[k],beta.prior.cov.sqrt[k,k]) ## always a 'valid' proposal
          beta.prop <- beta.curr
          beta.prop[k] <- bbeta
          
          
          mu.prop <- X %*% beta.prop
          
          llike.prop <- llike.GP(y, x, mu.prop, exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),2)
          
          lprior.prop <- logdmvnorm.sqrt(beta.prop, mu = beta.prior.mean, sqrt.Sigma = cov.beta.sqrt)
          
          l.acc <- llike.prop - llike.curr + lprior.prop - lprior.curr
          
          if(runif(1,0,1) <= min(l.acc,0)){
            beta.curr <- beta.prop
            
            llike.curr <- llike.prop
            
            lprior.curr <- lprior.prop
            
            mu.curr <- mu.prop
          }
          
        }else{
          ##use MH as normal
          
          ## first propose beta_k
          innovation <- rnorm(1,0,cov.beta.prop.sqrt[k,k])
          bbeta <- beta.curr[k] + innovation
          beta.prop <- beta.curr
          beta.prop[k] <- bbeta
          
          ## now accept/reject
          
          
          mu.prop <- X %*% beta.prop
          
          llike.prop <- llike.GP(y, x, mu.prop, exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),2)
          
          lprior.prop <- logdmvnorm.sqrt(beta.prop, mu = beta.prior.mean, sqrt.Sigma = cov.beta.sqrt)
          
          l.acc <- llike.prop - llike.curr + lprior.prop - lprior.curr
          
          if(log(runif(1,0,1)) <= min(l.acc,0)){
            beta.curr <- beta.prop
            
            llike.curr <- llike.prop
            
            lprior.curr <- lprior.prop
            
            mu.curr <- mu.prop
          }
          
          
          
          if(log(runif(1,0,1)) <= min(l.acc,0)){
            
            beta.curr <- beta.prop
            
            llike.curr <- llike.prop
            
            lprior.curr <- lprior.prop
            
            mu.curr <- mu.prop
          }
        }
        
      }
      
    }
    if(i > burn){
      ## store results
      
      ind <- i-burn
      
      psi.res[ind,] <- psi.curr
      beta.res[ind,] <- beta.curr
      llike.res[ind] <- llike.curr
      lprior.res[ind] <- lprior.curr
    }
  }
  list(psi = psi.res, beta = beta.res, llike = llike.res, lprior = lprior.res)
}


mcmc.ln.basis = function(y, x, basis.fns = id, nu=2, psi.curr, psi.mean, cov.psi.sqrt,
                         beta.curr, cov.beta.prop.sqrt, sqrt.prop.cov.mat,
                         beta.prior.mean, cov.beta.sqrt, p.prior, its=100, 
                         burn=0, thin=1){
  
  ## psi = (mu, log sig, log omega_vec, log lambda)
  beta.prior.cov.sqrt <- cov.beta.sqrt
  
  X <- basis.fns(x) ## this transforms x to h(x) (a vector of fns)
  
  n.psi <- length(psi.curr)
  n.beta <- length(beta.curr)
  mu.curr <- X %*% beta.curr #current value of the mean function at the design points
  
  
  llike.curr <- llike.GP(y,x,mu.curr,exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),nu)
  ## will assume log(hyper params) ~ N(.,.)
  lprior.curr <- logdmvnorm.sqrt(psi.curr,psi.mean,cov.psi.sqrt) + logdmvnorm.sqrt(beta.curr,beta.prior.mean,beta.prior.cov.sqrt)
  
  ## will work on log scale for hyper param
  
  logSig <- rep(0,its)
  n.Omega <- n.psi-2
  logOmega <- matrix(rep(0,its*n.Omega),ncol = n.Omega)
  logLambda <- rep(0,its)
  
  psi.res = cbind(logSig,logOmega,logLambda)
  
  ## beta will be on the standard scale
  
  n.beta = length(beta.curr)
  
  beta.res = matrix(rep(0,n.beta*its),ncol = n.beta, nrow = its)
  
  
  llike.res <- rep(0,its)
  lprior.res <- rep(0,its)
  
  vInnov.psi <- rep(0,n.psi)
  
  prop.cv.psi <- t(sqrt.prop.cov.mat) %*% sqrt.prop.cov.mat
  for(k in 1:n.psi){
    
    Sig12.p <- prop.cv.psi[-k,k]
    Sig21.p <- prop.cv.psi[k,-k]
    inv.sq <- solve(prop.cv.psi[-k,-k])
    
    tmp <- Sig12.p %*% inv.sq %*% Sig21.p
    #tmp2 <-  (inv.sq) %*% Sig21.p
    vInnov.psi[k] <- prop.cv.psi[k,k] - tmp #%*% tmp2
  }
  
  ## similar for beta
  
  vInnov.beta <- rep(0,n.beta)
  
  prop.cv.beta <- t(cov.beta.prop.sqrt) %*% cov.beta.prop.sqrt
  for(k in 1:n.beta){
    
    Sig12.p <- prop.cv.beta[-k,k]
    Sig21.p <- prop.cv.beta[k,-k]
    inv.sq <- solve(prop.cv.beta[-k,-k])
    
    tmp <- Sig12.p %*% inv.sq %*% Sig21.p
    #tmp2 <-  (inv.sq) %*% Sig21.p
    vInnov.beta[k] <- prop.cv.beta[k,k] - tmp #%*% tmp2
  }
  
  
  
  for(i in 1:(burn+its)){
    
    print(i-burn)
    
    for(j in 1:thin){
      
      ## first work through each component of psi
      for(k in 1:n.psi){
        ## propose from prior
        
        if(runif(1,0,1) <= p.prior){
          psi.prop <- psi.curr
          
          psi.prop[k] <- rnorm(1,psi.mean[k],cov.psi.sqrt[k,k])
          
          llike.prop <- llike.GP(y,x,mu.curr,exp(psi.prop[1]),exp(psi.prop[-c(1,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          lprior.prop <- logdmvnorm.sqrt(psi.prop,psi.mean,cov.psi.sqrt)
          
          l.acc <- llike.prop - llike.curr 
          
          if(log(runif(1,0,1))  <= min(l.acc,0)){
            psi.curr <- psi.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
          }
          
        }else{
          ## random walk prop
          
          psi.prop <- psi.curr
          
          psi.prop[k] <-    rnorm(1,psi.curr[k],sqrt(vInnov.psi[k]))
          
          
          ##NB this is always a 'valid' proposal
          
          
          llike.prop <- llike.GP(y,x,mu.curr,exp(psi.prop[1]),exp(psi.prop[-c(1,n.psi)]),exp(psi.prop[n.psi]),nu)
          
          lprior.prop <- logdmvnorm.sqrt(psi.prop, mu = psi.mean, sqrt.Sigma = cov.psi.sqrt)
          ## now eval acc prob
          
          l.acc <- llike.prop - lprior.curr + lprior.prop - llike.curr ## sym chain therefore do not consider q(.,.)
          
          if( log(runif(1,0,1)) <= min(l.acc,0)){
            psi.curr <- psi.prop
            llike.curr <- llike.prop
            lprior.curr <- lprior.prop
          }
          
        }
        
      }
      ##now for the betas which will again be done in steps
      for(k in 1:n.beta){
        ## first propose from prior
        if(runif(1,0,1) <= p.prior){
          
          ## prop the k-th component of beta
          
          bbeta <- rnorm(1,beta.prior.mean[k],beta.prior.cov.sqrt[k,k]) ## always a 'valid' proposal
          
          ## NB in prior beta~iid a prior - no need to condition
          beta.prop <- beta.curr
          beta.prop[k] <- bbeta
          
          
          mu.prop <- X %*% beta.prop
          
          llike.prop <- llike.GP(y, x, mu.prop, exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),2)
          
          lprior.prop <- logdmvnorm.sqrt(beta.prop, mu = beta.prior.mean, sqrt.Sigma = cov.beta.sqrt)
          
          l.acc <- llike.prop - llike.curr + lprior.prop - lprior.curr
          
          if(log(runif(1,0,1)) <= min(l.acc,0)){
            beta.curr <- beta.prop
            
            llike.curr <- llike.prop
            
            lprior.curr <- lprior.prop
            
            mu.curr <- mu.prop
          }
          
        }else{
          ##use MH as normal
          
          ## first propose beta_k
          
          ## seperate out the covariance
          
          innovation <- rnorm(1,0,vInnov.beta^0.5)
          bbeta <- beta.curr[k] + innovation
          beta.prop <- beta.curr
          beta.prop[k] <- bbeta
          
          ## now accept/reject
          
          
          mu.prop <- X %*% beta.prop
          
          llike.prop <- llike.GP(y, x, mu.prop, exp(psi.curr[1]),exp(psi.curr[-c(1,n.psi)]),exp(psi.curr[n.psi]),2)
          
          lprior.prop <- logdmvnorm.sqrt(beta.prop, mu = beta.prior.mean, sqrt.Sigma = cov.beta.sqrt)
          
          l.acc <- llike.prop - llike.curr + lprior.prop - lprior.curr
          
          if(log(runif(1,0,1)) <= min(l.acc,0)){
            beta.curr <- beta.prop
            
            llike.curr <- llike.prop
            
            lprior.curr <- lprior.prop
            
            mu.curr <- mu.prop
          }
          
          
          
          #if(log(runif(1,0,1)) <= min(l.acc,0)){
          
          # beta.curr <- beta.prop
          
          #llike.curr <- llike.prop
          
          #lprior.curr <- lprior.prop
          
          #mu.curr <- mu.prop
          #}
        }
        
      }
      
    }
    if(i > burn){
      ## store results
      
      ind <- i-burn
      
      psi.res[ind,] <- psi.curr
      beta.res[ind,] <- beta.curr
      llike.res[ind] <- llike.curr
      lprior.res[ind] <- lprior.curr
    }
  }
  list(psi = psi.res, beta = beta.res, llike = llike.res, lprior = lprior.res)
}
