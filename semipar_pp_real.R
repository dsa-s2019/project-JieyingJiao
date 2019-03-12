rm(list=ls())
library(nimble)
library(ggplot2)
library(rjson)
library(grid)
library(jpeg)
library(RCurl)
library(coda)
library(gridExtra)
library(reshape2)
library(matrixStats)
library(viridis)

###=============================================================================

### data visualization

##=============================================================================


## analysis real data, cut half court, only focuse on the field that have most shots
## find the max y lab of a made shot
name.vector <- list.files("real_data/")
shotdata <- vector("list", length = length(name.vector))
y.max <- rep(0, length = length(name.vector))
for (i in 1:length(name.vector)) {
  shotdata[[i]] <- read.csv(paste0("real_data/", name.vector[i]))
  shotdata[[i]] <- shotdata[[i]][as.logical((shotdata[[i]]$y <= 420)*(shotdata[[i]]$y >= -7.5)), ]
  y.max[i] <- max(shotdata[[i]]$y[shotdata[[i]]$shot_made_flag == 1])
}

## most of the made shot are below y = 300, only a few are above that
## so cut the half court at y = 300
## original scale is [-250, 250] * [0, 300]
for (i in 1:length(name.vector)) {
  shotdata[[i]] <- shotdata[[i]][shotdata[[i]]$y < 300, ]
}


## choose one person's data to fit
name <- 5
name.vector[name]

ggplot(data = shotdata[[name]], aes(x = x, y = y)) + 
  geom_point(aes(colour = as.factor(shot_made_flag)))


## full model
## intensity model: x, y, distance
## convert x, y into [0, 25] * [0, 30] scale
## symmetric in x about basket rim
x <- shotdata[[name]]$x / 10
y <- shotdata[[name]]$y / 10

## creat dummy variables, 0 for 2PT, 1 for 3PT
shot_type <- model.matrix(~shotdata[[name]]$shot_type)[, -1]
period <- model.matrix(~factor(shotdata[[name]]$period))[, -1]

## marks model: intensity, intercept, period(treat as integer), 
##              seconds_left(minutes*60+seconds), shot_type(factor), distance
seconds <- shotdata[[name]]$minutes_remaining * 60 + shotdata[[name]]$seconds_remaining


##' define log-likelihood
llFun <- nimbleFunction(
  setup <- function(model) {},
  
  run <- function() {
    s_ll <- model$s_ll
    log_ll <- model$log_ll
    ll <- log_ll - s_ll
    returnType(double())
    return(ll[1])
  }
)



## shot type index matrix (0 for 2PT) and d.adv matrix on [0, 1] square window 
## according to only right half
index.type <- matrix(1, nrow = 123, ncol = 200)
dadv <- matrix(0, nrow = 123, ncol = 200)
angle1 <- matrix(0, nrow = 123, ncol = 200)
angle2 <- matrix(0, nrow = 123, ncol = 200)
angle3 <- matrix(0, nrow = 123, ncol = 200)
angle4 <- matrix(0, nrow = 123, ncol = 200)
angle5 <- matrix(0, nrow = 123, ncol = 200)
for (i in 1:200) {
  for (j in 1:123) {
    
    x.o <- (i-100)/4 * 10
    y.o <- (j-3)/4 * 10
    d.o <- floor(sqrt(x.o^2 + y.o^2)/10)
    
    if ((y.o<=90)*(abs(x.o)<=220)) {
      index.type[j, i] <- 0
    } else if((y.o>90)*(y.o<=237.5)*(sqrt(x.o^2+y.o^2)<=237.5)) index.type[j, i] <- 0
    
    dadv[j, i] <- (1 - index.type[j, i]) * d.o + 
      index.type[j, i] * (y.o<=90) * floor((abs(x.o)-220)/10) + 
      index.type[j, i] * (y.o>90) * floor((sqrt(x.o^2+y.o^2) - 237.5)/10)
    
    angle.value <- atan2(y.o, x.o)
    if (angle.value < -pi/2) angle.value <- 2*pi+angle.value
    if (angle.value <= pi/6) {
      angle1[j, i] <- 0
      angle2[j, i] <- 0
      angle3[j, i] <- 0
      angle4[j, i] <- 0
      angle5[j, i] <- 0
    } else if (angle.value <= pi/3) {
      angle1[j, i] <- 1
      angle2[j, i] <- 0
      angle3[j, i] <- 0
      angle4[j, i] <- 0
      angle5[j, i] <- 0
    } else if (angle.value <= pi/2) {
      angle1[j, i] <- 0
      angle2[j, i] <- 1
      angle3[j, i] <- 0
      angle4[j, i] <- 0
      angle5[j, i] <- 0
    } else if (angle.value <= 2*pi/3) {
      angle1[j, i] <- 0
      angle2[j, i] <- 0
      angle3[j, i] <- 1
      angle4[j, i] <- 0
      angle5[j, i] <- 0
    } else if (angle.value <= 5*pi/6) {
      angle1[j, i] <- 0
      angle2[j, i] <- 0
      angle3[j, i] <- 0
      angle4[j, i] <- 1
      angle5[j, i] <- 0
    } else {
      angle1[j, i] <- 0
      angle2[j, i] <- 0
      angle3[j, i] <- 0
      angle4[j, i] <- 0
      angle5[j, i] <- 1
    }
    
  }
}



## standardize dadv
dadv <- (dadv-mean(dadv)) / sd(dadv)


constants <- list(N = nrow(dadv)*ncol(dadv), M = 20, 
                  Nshot = nrow(shotdata[[name]]))

data <- list(
  index.type = index.type,
  dadv = dadv, 
  rowindex = floor(4*(y+0.75)) + 1,
  colindex = floor(4*(x+25)) + 1
)


###============================================================================

### variable selection using spike-slab prior
### fit and select intensity dependent model

###============================================================================

code <- nimbleCode({
  beta1 ~ dnorm(0, sd = 10) # shot_type
  beta2 ~ dnorm(0, sd = 10) # dadv
  beta3 ~ dnorm(0, sd = 10) # shot_type * dadv
  
  
  z[1:N] ~ dCRP(alpha, size = N)
  alpha ~ dgamma(1, 1)
  for (i in 1:M) {
    lambda0[i] ~ dgamma(1, 1)
  }
  

  for (i in 1:200) {
    for (j in 1:123) {
      lambda_D[j, i] <- lambda0[z[(j-1)*200+i]] * 
        exp(beta1 * index.type[j, i] + beta2 * dadv[j, i] + 
              beta3 * index.type[j, i] * dadv[j, i])
    }
  }
  
  
  s_ll <- sum(lambda_D[1:123, 1:200])/16
  
  for (i in 1:Nshot) {
    lambda[i] <- lambda0[z[(rowindex[i]-1)*200+colindex[i]]] * 
      exp(beta1 * index.type[rowindex[i], colindex[i]] + 
            beta2 * dadv[rowindex[i], colindex[i]] + 
            beta3 * index.type[rowindex[i], colindex[i]] * 
            dadv[rowindex[i], colindex[i]])
    
  }
  log_ll <- sum(log(lambda[1:Nshot]))
})

## initial value of MCMC sampling
inits <- list(beta1 = 0, beta2 = 0, beta3 = 0, lambda0 = rep(1, constants$M), 
              z = sample(1:10, size = constants$N, replace = TRUE), 
              alpha  = 1)

Rmodel <- nimbleModel(code = code, constants = constants, data = data, 
                      inits = inits, check = FALSE)
RllFun <- llFun(Rmodel)
mcmcConf <- configureMCMC(Rmodel)

## loglikelihood control MH algorithm
mcmcConf$addMonitors(target = c(paste0("beta", 1:3), "lambda0", "z", "alpha"))

for (i in 1:3) {
  name <- paste0("beta", i)
  mcmcConf$addSampler(target = name, type = 'RW_llFunction',
                      control = list(llFunction = RllFun,
                                     includesTarget = FALSE)) 
}


for (i in constants$M) {
  mcmcConf$addSampler(target = paste0("lambda0[", i, "]"), type = 'RW_llFunction',
                      control = list(llFunction = RllFun,
                                     includesTarget = FALSE))
}

## Complie Command
Rmcmc <- buildMCMC(mcmcConf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
## MCMC run
Cmcmc$run(60000)
## get MCMC sample
samples <- as.matrix(Cmcmc$mvSamples)[seq(20001, 60000, by = 10), ]

plot(as.mcmc(samples))