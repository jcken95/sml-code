library("rstan")
source("../GPfunctionsOptim.R")
source("../hetGPfunctions.R")
het.uni <- stan_model("../hetGP-univariate.stan")
SML <- stan_model("../SML-hetGP-univariate.stan")
## lets define some functions ...

f1 <- function(x){
  4*sin(6*pi*x)
}
f2 <- function(x){
  4*sin(6*pi*x) + 5*(2*x+1)
}
f1.rand <- function(x){
  
  f1(x) + 0.5*rnorm(length(x))
  
}

f2.rand <- function(x){
  
  f2(x) + (1.1 + 7*x)*rnorm(length(x))
  
}
plot(x, f2(x))
plot(x,f2.rand(x))
plot(f1(x), f2(x))
NN <- 100 ## number reps for sim study

x <- seq(0 , 1, length = NN)
plot(x, f1(x), type = "l")
plot(x, f2.rand(x), type = "l")
plot(f1(x), f2(x)); cor(f1(x), f2(x))



############################
############################
############################
######## start here ########
############################
############################
############################
NN <- 100
score.emp1 <- mse.emp1 <- matrix(0, ncol = 2, nrow = NN)
set.seed(1234)
st1 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = rep(0, 2), m_beta_s = rep(10, 2),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = rep(0, 2), v_beta_s = rep(10, 2),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  Var.x <- cond.var.het(pars$logLambda, x.std, x.new, pars$v_sigma,
                        pars$v_theta, pars$v_nugget, as.vector(cbind(1, x.new)%*%pars$v_beta), 
                        as.vector(cbind(1, x.std)%*%pars$v_beta), 2)
  Var.x <- diag( cond.covmat.het(exp(pars$logLambda), Var.x, x.std, x.new, pars$m_sigma, pars$m_theta, 2) )
  
  Mean.x <- cond.mean.het(y1, x.std, x.new, pars$m_sigma, pars$m_theta, exp(pars$logLambda), 1e-8,
                          cbind(1, x.new)%*%pars$m_beta, cbind(1, x.std)%*%pars$m_beta, 2)
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  n.c <- 100
  n.e <- length(x1) - 1
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = rep(0, 2), m_beta_s = rep(10,2),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = rep(0, 2), v_beta_s = rep(10,2),
    v_a_theta = 2, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = rep(0, 2), c_beta_s = rep(10,2),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 1, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  Var <- rep(0, length(x.new))
  Var <- cond.var.hetSML(pars$logLambda, x2, x.new, 
                         pars$v_sigma, pars$v_theta,  pars$v_nugget, 
                         as.vector(cbind(1, x.new)%*%pars$v_beta), 
                         as.vector(cbind(1, x2)%*%pars$v_beta),2)
  Var
  
  Var <- diag( cond.covmat.hetSML(exp(pars$logLambda), Var,  c(x.c, x2), x.new,
                                  pars$m_sigma, pars$m_theta, pars$c_sigma, pars$c_theta, pars$c_nugget,
                                  pars$rho, n.c, n.c + n.e, 2)  )
  mu.yobs <- as.vector(rbind(cbind(1, x.c)%*%pars$c_beta, cbind(1, x2)%*%(pars$rho*pars$c_beta + pars$m_beta)))
  mu.y <- as.vector(cbind(1, x.new)%*%(pars$rho*pars$c_beta + pars$m_beta))
  Mean <- cond.mean.hetSML(c(y.c, y2), c(x.c, x2), x.new, c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                           c(pars$m_sigma, pars$m_theta, 0),
                           exp(pars$logLambda), pars$rho, mu.y, mu.yobs, n.c, 2)
  test.data <- f2.rand(x.test)
  
  mse.emp1[K, ] <-   c(  mean((Mean - test.data)^2), mean((Mean.x - test.data)^2) )
  score.emp1[K, ] <- c(sum( (-(Mean - test.data)^2)/Var - log(Var) ), sum( (-(Mean.x - test.data)^2)/Var.x - log(Var.x)))
  
  
  print(K)
  
}
en1 <- Sys.time()
tot1 <- en1 - st1

tot1
boxplot(mse.emp1)
boxplot(score.emp1, ylim = c(-4400,-3800))
#### simulation 2 ####


score.emp2 <- mse.emp2 <- matrix(0, ncol = 2, nrow = NN)

st2 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = rep(0, 2), m_beta_s = rep(10, 2),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = rep(0, 2), v_beta_s = rep(10, 2),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  Var.x <- cond.var.het(pars$logLambda, x.std, x.new, pars$v_sigma,
                        pars$v_theta, pars$v_nugget, as.vector(cbind(1, x.new)%*%pars$v_beta), 
                        as.vector(cbind(1, x.std)%*%pars$v_beta), 2)
  Var.x <- diag( cond.covmat.het(exp(pars$logLambda), Var.x, x.std, x.new, pars$m_sigma, pars$m_theta, 2) )
  
  Mean.x <- cond.mean.het(y1, x.std, x.new, pars$m_sigma, pars$m_theta, exp(pars$logLambda), 1e-8,
                          cbind(1, x.new)%*%pars$m_beta, cbind(1, x.std)%*%pars$m_beta, 2)
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  n.c <- 200
  n.e <- length(x1) - 2
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = rep(0, 2), m_beta_s = rep(10,2),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = rep(0, 2), v_beta_s = rep(10,2),
    v_a_theta = 5, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = rep(0, 2), c_beta_s = rep(10,2),
    c_a_theta = 5, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 1, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  Var <- rep(0, length(x.new))
  Var <- cond.var.hetSML(pars$logLambda, x2, x.new, 
                         pars$v_sigma, pars$v_theta,  pars$v_nugget, 
                         as.vector(cbind(1, x.new)%*%pars$v_beta), 
                         as.vector(cbind(1, x2)%*%pars$v_beta),2)
  Var
  
  Var <- diag( cond.covmat.hetSML(exp(pars$logLambda), Var,  c(x.c, x2), x.new,
                                  pars$m_sigma, pars$m_theta, pars$c_sigma, pars$c_theta, pars$c_nugget,
                                  pars$rho, n.c, n.c + n.e, 2)  )
  mu.yobs <- as.vector(rbind(cbind(1, x.c)%*%pars$c_beta, cbind(1, x2)%*%(pars$rho*pars$c_beta + pars$m_beta)))
  mu.y <- as.vector(cbind(1, x.new)%*%(pars$rho*pars$c_beta + pars$m_beta))
  Mean <- cond.mean.hetSML(c(y.c, y2), c(x.c, x2), x.new, c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                           c(pars$m_sigma, pars$m_theta, 0),
                           exp(pars$logLambda), pars$rho, mu.y, mu.yobs, n.c, 2)
  test.data <- f2.rand(x.test)
  
  mse.emp2[K, ] <-   c(  mean((Mean - test.data)^2), mean((Mean.x - test.data)^2) )
  score.emp2[K, ] <- c(sum( (-(Mean - test.data)^2)/Var - log(Var) ), sum( (-(Mean.x - test.data)^2)/Var.x - log(Var.x)))
  
  
  print(K)
  
}
en2 <- Sys.time()
tot2 <- en2 - st2

tot2
boxplot(mse.emp2, ylim=c(9,15))
boxplot(score.emp2, ylim= c(-4000,-3000))

boxplot(mse.emp1, ylim=c(9,15))
boxplot(score.emp1, ylim= c(-4000,-3000))

score.emp3 <- mse.emp3 <- matrix(0, ncol = 2, nrow = NN)

st3 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = rep(0, 2), m_beta_s = rep(10, 2),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = rep(0, 2), v_beta_s = rep(10, 2),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  Var.x <- cond.var.het(pars$logLambda, x.std, x.new, pars$v_sigma,
                        pars$v_theta, pars$v_nugget, as.vector(cbind(1, x.new)%*%pars$v_beta), 
                        as.vector(cbind(1, x.std)%*%pars$v_beta), 2)
  Var.x <- diag( cond.covmat.het(exp(pars$logLambda), Var.x, x.std, x.new, pars$m_sigma, pars$m_theta, 2) )
  
  Mean.x <- cond.mean.het(y1, x.std, x.new, pars$m_sigma, pars$m_theta, exp(pars$logLambda), 1e-8,
                          cbind(1, x.new)%*%pars$m_beta, cbind(1, x.std)%*%pars$m_beta, 2)
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  n.c <- 300
  n.e <- length(x1) - 3
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = rep(0, 2), m_beta_s = rep(10,2),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = rep(0, 2), v_beta_s = rep(10,2),
    v_a_theta = 2, v_b_theta = 1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = rep(0, 2), c_beta_s = rep(10,2),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 1, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  Var <- rep(0, length(x.new))
  Var <- cond.var.hetSML(pars$logLambda, x2, x.new, 
                         pars$v_sigma, pars$v_theta,  pars$v_nugget, 
                         as.vector(cbind(1, x.new)%*%pars$v_beta), 
                         as.vector(cbind(1, x2)%*%pars$v_beta),2)
  Var
  
  Var <- diag( cond.covmat.hetSML(exp(pars$logLambda), Var,  c(x.c, x2), x.new,
                                  pars$m_sigma, pars$m_theta, pars$c_sigma, pars$c_theta, pars$c_nugget,
                                  pars$rho, n.c, n.c + n.e, 2)  )
  mu.yobs <- as.vector(rbind(cbind(1, x.c)%*%pars$c_beta, cbind(1, x2)%*%(pars$rho*pars$c_beta + pars$m_beta)))
  mu.y <- as.vector(cbind(1, x.new)%*%(pars$rho*pars$c_beta + pars$m_beta))
  Mean <- cond.mean.hetSML(c(y.c, y2), c(x.c, x2), x.new, c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                           c(pars$m_sigma, pars$m_theta, 0),
                           exp(pars$logLambda), pars$rho, mu.y, mu.yobs, n.c, 2)
  test.data <- f2.rand(x.test)
  
  mse.emp3[K, ] <-   c(  mean((Mean - test.data)^2), mean((Mean.x - test.data)^2) )
  score.emp3[K, ] <- c(sum( (-(Mean - test.data)^2)/Var - log(Var) ), sum( (-(Mean.x - test.data)^2)/Var.x - log(Var.x)))
  
  
  print(K)
  
}
en3 <- Sys.time()
tot3 <- en3 - st3



score.emp4 <- mse.emp4 <- matrix(0, ncol = 2, nrow = NN)

st4 <- Sys.time()
for(K in 1:NN){
  
  ## first generate "purely exp" emulator
  
  
  x1 <- lhs::maximinLHS(50, 1)
  #plot(x1, f2.rand(x1))
  lines(x, f2(x))
  #x1 <- lhs::maximinLHS(10, 1)
  y1 <- f2.rand(x1)
  #plot(x1, y1)
  ## test x
  x.test <- as.vector(lhs::maximinLHS(1000,1))
  x.std <- scale(x1) 
  dist.x <- as.matrix(dist(x.std, upper = T, diag = T))
  data.het <- list(
    m_p = 2, v_p = 2,  N = length(x1), K = 1,
    x = as.matrix(x.std), y = as.vector(y1), Diffs = dist.x,
    m_H = cbind(1, x.std), v_H = cbind(1, x.std),
    ## prior ##
    
    m_beta_m = rep(0, 2), m_beta_s = rep(10, 2),
    
    m_a_theta = 2, m_b_theta = 1,
    
    m_a_sigma = 2, m_b_sigma = 2,
    
    m_nugget = 0,
    
    v_beta_m = rep(0, 2), v_beta_s = rep(10, 2),
    
    v_a_theta = 2, v_b_theta = 1,
    
    v_a_sigma = 2, v_b_sigma = 2,
    
    v_nugget_a = 2, v_nugget_b = 2
  )
  
  het.uni.fit <- rstan::optimizing(het.uni, data.het, as_vector = F)
  het.uni.fit
  
  pars <- het.uni.fit$par
  
  ## plot resulting emulator ...
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  Var.x <- cond.var.het(pars$logLambda, x.std, x.new, pars$v_sigma,
                        pars$v_theta, pars$v_nugget, as.vector(cbind(1, x.new)%*%pars$v_beta), 
                        as.vector(cbind(1, x.std)%*%pars$v_beta), 2)
  Var.x <- diag( cond.covmat.het(exp(pars$logLambda), Var.x, x.std, x.new, pars$m_sigma, pars$m_theta, 2) )
  
  Mean.x <- cond.mean.het(y1, x.std, x.new, pars$m_sigma, pars$m_theta, exp(pars$logLambda), 1e-8,
                          cbind(1, x.new)%*%pars$m_beta, cbind(1, x.std)%*%pars$m_beta, 2)
  
  
  ## can we add in some points from sin(4pix) to help?
  
  
  ## can we build a multilevel emulator instead??
  n.c <- 400
  n.e <- length(x1) - 4
  x2 <- lhs::maximinLHS(n.e,1)
  y2 <- f2.rand(x2)
  x.c <- rev(lhs::augmentLHS(x2, n.c-n.e))
  x.c.old <- x.c
  x2 <- rev(x2)
  
  y.c <- f1.rand(x.c)
  x.c <- scale(x.c, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  x2 <- scale(x2, center = attr(x.std,"scaled:center"), scale = attr(x.std,"scaled:scale"))
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  dists <- as.matrix(dist(x.c, diag = T, upper = T))
  y.c <- y.c; y2 <- rev(y2)
  H <- matrix(0, ncol = 4, nrow = n.c)
  H[,1:2] <- cbind(1, x.c)
  H[(n.c - n.e + 1):n.c, 3:4] <- cbind(1, x2)
  
  #x.c <- scale(x.c, center = attr(x2,"scaled:center"), scale = attr(x2,"scaled:scale"))
  #dists <- as.matrix(dist(x.c, diag = T, upper = T))
  sml.data <- list(
    m_p = 2, v_p = 2, N = n.e + n.c,
    n_c = n.c, n_e = n.e,
    x_e = as.matrix(x2), x_c = as.matrix(x.c),
    v_H = cbind(1, x2), m_H = H,
    y_c = as.vector(y.c), y_e = as.vector(y2),
    ## now prior
    
    # mean
    m_beta_m = rep(0, 2), m_beta_s = rep(10,2),
    m_a_theta = 2, m_b_theta = 1, 
    m_a_sigma = 2, m_b_sigma = 2,
    m_nugget = 0,
    
    # variance
    v_beta_m = rep(0, 2), v_beta_s = rep(10,2),
    v_a_theta = 2, v_b_theta =  1, 
    v_a_sigma = 2, v_b_sigma = 2,
    v_nugget_a = 2, v_nugget_b = 2,
    
    # cheap code
    c_beta_m = rep(0, 2), c_beta_s = rep(10,2),
    c_a_theta = 2, c_b_theta = 1, 
    c_a_sigma = 2, c_b_sigma = 2,
    c_nugget_a = 2, c_nugget_b = 2,
    
    D_cc = dists, D_ee = dists[(n.c - n.e + 1):n.c, (n.c - n.e + 1):n.c],
    D_ce = dists[, (n.c - n.e + 1):n.c],
    
    m_rho = 1, s_rho = 1/3
    
  )
  sml.fit <- rstan::optimizing(SML, sml.data, as_vector = F)
  sml.fit
  
  pars <- sml.fit$par
  x.new <- scale(x.test, center = attr(x.std,"scaled:center"), scale =  attr(x.std,"scaled:scale"))
  Var <- rep(0, length(x.new))
  Var <- cond.var.hetSML(pars$logLambda, x2, x.new, 
                         pars$v_sigma, pars$v_theta,  pars$v_nugget, 
                         as.vector(cbind(1, x.new)%*%pars$v_beta), 
                         as.vector(cbind(1, x2)%*%pars$v_beta),2)
  Var
  
  Var <- diag( cond.covmat.hetSML(exp(pars$logLambda), Var,  c(x.c, x2), x.new,
                                  pars$m_sigma, pars$m_theta, pars$c_sigma, pars$c_theta, pars$c_nugget,
                                  pars$rho, n.c, n.c + n.e, 2)  )
  mu.yobs <- as.vector(rbind(cbind(1, x.c)%*%pars$c_beta, cbind(1, x2)%*%(pars$rho*pars$c_beta + pars$m_beta)))
  mu.y <- as.vector(cbind(1, x.new)%*%(pars$rho*pars$c_beta + pars$m_beta))
  Mean <- cond.mean.hetSML(c(y.c, y2), c(x.c, x2), x.new, c(pars$c_sigma, pars$c_theta, pars$c_nugget),
                           c(pars$m_sigma, pars$m_theta, 0),
                           exp(pars$logLambda), pars$rho, mu.y, mu.yobs, n.c, 2)
  test.data <- f2.rand(x.test)
  
  mse.emp4[K, ] <-   c(  mean((Mean - test.data)^2), mean((Mean.x - test.data)^2) )
  score.emp4[K, ] <- c(sum( (-(Mean - test.data)^2)/Var - log(Var) ), sum( (-(Mean.x - test.data)^2)/Var.x - log(Var.x)))
  
  
  print(K)
  
}
en4 <- Sys.time()
tot4 <- en4 - st4
boxplot(mse.emp4)
boxplot(cbind(score.emp1, score.emp2, score.emp3, score.emp4), ylim = c(-5000, -3800))
boxplot(cbind(mse.emp1, mse.emp2, mse.emp3, mse.emp4))
