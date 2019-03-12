rm(list=ls())
library(nimble)
library(spatstat)


##' data simulation
##' lamnbda0 is a vector for different clusters, now focuse on 4 clusters
##' beta is a two dimension vector for (x, y) coefficients
##' grid = c(nrow, ncol)
##' win1, win2 are two dimension vectors showing area to generate data
data_simu <- function(lambda0, beta, grid, win1, win2) {
  x.grid <- seq(win1[1], by = (win1[2]-win1[1])/grid[2], length = grid[2])
  y.grid <- seq(win2[1], by = (win2[2] - win2[1])/grid[1], length = grid[1])
  
  lambda.grid <- matrix(0, nrow = grid[1], ncol = grid[2])
  
  for (i in 1:floor(grid[1]/2)) {
    for (j in 1:floor(grid[2]/2)) {
      lambda.grid[i, j] <- lambda0[1] * 
        exp(beta[1] * x.grid[j] + beta[2] * y.grid[i])
    }
    for (j in (floor(grid[2]/2)+1):grid[2]) {
      lambda.grid[i, j] <- lambda0[2] * 
        exp(beta[1] * x.grid[j] + beta[2] * y.grid[i])
    }
  }
  
  for (i in (floor(grid[1]/2)+1) : grid[1]) {
    for (j in 1:floor(grid[2]/2)) {
      lambda.grid[i, j] <- lambda0[3] * 
        exp(beta[1] * x.grid[j] + beta[2] * y.grid[i])
    }
    for (j in (floor(grid[2]/2)+1):grid[2]) {
      lambda.grid[i, j] <- lambda0[4] * 
        exp(beta[1] * x.grid[j] + beta[2] * y.grid[i])
    }
  }
  
  pp <- rpoispp(as.im(lambda.grid, W = owin(win1, win2)))
  
  N.grid <- matrix(0, nrow = grid[1], ncol = grid[2])
  x.grid <- c(x.grid, win1[2])
  y.grid <- c(y.grid, win2[2])
  for (i in 1:grid[1]) {
    for (j in 1:grid[2]) {
      N.grid[i, j] <- sum((pp$x>x.grid[j]) * (pp$x<x.grid[j+1]) * 
                      (pp$y>y.grid[i]) * (pp$y<y.grid[i+1]))
    }
  }
  
  return(list(x = pp$x, y = pp$y, N.grid = N.grid))
}

##' define log-likelihood
llFun <- nimbleFunction(
  setup <- function(model) {},
  
  run <- function() {
    ll <- model$log_ll
    returnType(double())
    return(ll[1])
  }
)

lambda0 <- c(10, 10, 5000, 1000)
beta <- c(0, 0)
grid <- c(20, 20)
win1 <- c(-1, 1)
win2 <- c(-1, 1)


simudata <- data_simu(lambda0 = lambda0, beta = beta, grid = grid, 
                      win1 = win1, win2 = win2)

constants <- list(grid = c(20, 20), N = 400, M = 20)


data <- list(
  N.grid = simudata$N.grid,
  x.grid = seq(win1[1], by = (win1[2]-win1[1])/grid[2], length = grid[2]),
  y.grid = seq(win2[1], by = (win2[2] - win2[1])/grid[1], length = grid[1])
)



code <- nimbleCode({
  
  
  z[1:N] ~ dCRP(alpha, size = N)
  alpha ~ dgamma(1, 1)
  for (i in 1:M) {
    lambda0[i] ~ dgamma(1, 1)
  }
  
  for (i in 1:grid[2]) {
    for (j in 1:grid[1]) {
      lambda[j, i] <- lambda0[z[(j-1)*grid[2]+i]] 
      N.grid[j, i] ~ dpois(lambda[j, i])
    }
  }
  
})

inits <- list(lambda0 = rep(1, constants$M), 
              z = sample(1:10, size = 400, replace = TRUE), 
              alpha  = 1)

Rmodel <- nimbleModel(code = code, constants = constants, data = data, 
                      inits = inits, check = FALSE)
#RllFun <- llFun(Rmodel)
mcmcConf <- configureMCMC(Rmodel)

mcmcConf$addMonitors(target = c( "lambda0", "z", "alpha"))


Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
## MCMC run
Cmcmc$run(10000)
## get MCMC sample
samples <- as.matrix(Cmcmc$mvSamples)[5001:10000, ]

plot(as.mcmc(samples))
