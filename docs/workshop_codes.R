

##
## Session 1: Introduction to NIMBLE
##

library(nimble)

code <- nimbleCode({
    for(i in 1:G) {
        ## priors:
        a[i] ~ dgamma(1, 0.001)
        b[i] ~ dgamma(1, 0.001)
        ##
        for(j in 1:N) {
            ## random effects survival
            p[i,j] ~ dbeta(a[i], b[i])
            ## likelihood:
            r[i,j] ~ dbinom(prob = p[i,j],
                            size = n[i,j])
        }
    }
})


load('~/Downloads/litters_data.Rdata')

G <- dim(n)[1]
N <- dim(n)[2]

constants <- list(N = N, G = G, n = n)

data <- list(r = r)

inits <- list(a = c(1, 1),
              b = c(1, 1),
              p = array(0.5, c(G,N)))


Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$a <- c(10, -10)
Rmodel$a

Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)

set.seed(0)
Rmodel$a
Rmodel$simulate('a')
calculate(Rmodel)

Cmodel$a

calculate(Rmodel)
calculate(Cmodel)

Rmodel$a <- c(1, 2)
Rmodel$a
Rmodel$logProb_a
Rmodel$calculate('a')
Rmodel$getLogProb('a')

Rmodel$getNodeNames(stochOnly = TRUE, includeData = FALSE)

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$calculate()  ##  -173.7622

conf <- configureMCMC(Rmodel)

## monitors
conf$printMonitors()
conf$setThin(10)
conf$setThin(1)
conf$printMonitors()
conf$addMonitors('p')

## samplers
conf$printSamplers()

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)



## let's you specify
## burnin, thinning, intial values,
## number of chains, progress bar
## samples as a coda object
## only a summary
samples <- runMCMC(Rmcmc, 200)
samples

set.seed(0)
samples <- runMCMC(Cmcmc, 20000)
samples

dim(samples)
colnames(samples)

## provided with nimble
round(samplesSummary(samples), 3)

samples <- runMCMC(Cmcmc, 20000, nburnin = 10000)

dim(samples)
str(out)


samplesList <- runMCMC(Cmcmc, 50000, nchains = 3)
str(samplesList)

initsFunction <- function() {
    list(a = runif(2, 0, 100),
         b = runif(2, 0, 100))
}

samplesList <- runMCMC(Cmcmc, 50000, nchains=3,
                       inits = initsFunction)


## effective sample size:

library(coda)

effectiveSize(rnorm(10000))
effectiveSize(c(1:1000))
apply(samples, 2, effectiveSize)

samples <- samplesList[[1]]


## library basicMCMCplots:

library(basicMCMCplots)

samplesPlot(samples, 'a')
samplesPlot(samples, 'b')

chainsPlot(samples)
chainsPlot(samplesList)



##
## Session 2: Improving MCMC Performance, and Dipper Model
##


code <- nimbleCode({
    for (i in 1:G) {
        a[i] ~ dgamma(1, .001)
        b[i] ~ dgamma(1, .001)
        for (j in 1:N) {
            ## random effects
            p[i,j] ~ dbeta(a[i], b[i]) 
            ## likelihood
            r[i,j] ~ dbin(p[i,j], n[i,j])
        }
    }
})


## slice sampling for litters model:

load('~/Downloads/litters_data.Rdata')

G <- dim(n)[1]
N <- dim(n)[2]

constants <- list(N = N, G = G, n = n)

data <- list(r = r)

inits <- list(a = c(1, 1),
              b = c(1, 1),
              p = array(0.5, c(G,N)))

Rmodel <- nimbleModel(code, constants,
                      data, inits)

conf <- configureMCMC(Rmodel)
conf <- configureMCMC(Rmodel, onlySlice = TRUE)

conf$printSamplers()

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 10000)

round(samplesSummary(samples), 2)

dim(samples)
colnames(samples)
samplesPlot(samples[1:1000, 1:2], densityplot = FALSE)


library(coda)

apply(samples, 2, effectiveSize)

## what is happening here?
## basically, it's due to
## high cross-correlation in the a's and b's

library(basicMCMCplots)

colnames(samples)

samplesPlot(samples[, c(2,4)], scale = TRUE, burnin = 8000)

round(cor(samples), 2)


## block sampling for litters model:

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)
conf$printSamplers()

conf$removeSamplers(c('a', 'b'))
conf$printSamplers()

conf$addSampler(target = c('a[1]', 'b[1]'),
                type = 'RW_block')

conf$printSamplers()

conf$addSampler(target = c('a[2]', 'b[2]'),
                type = 'RW_block')

conf$printSamplers()

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

system.time(samples <- runMCMC(Cmcmc, 10000))

cor(samples)
cov(samples)

apply(samples, 2, effectiveSize)
round(samplesSummary(samples), 2)


## block sampling with initial covariance matrices:

Rmcmc$samplerFunctions$contentsList[[33]]
Rmcmc$samplerFunctions$contentsList[[33]]$target
Rmcmc$samplerFunctions$contentsList[[33]]$adaptive
Rmcmc$samplerFunctions$contentsList[[33]]$adaptInterval
Rmcmc$samplerFunctions$contentsList[[33]]$propCov
Rmcmc$samplerFunctions$contentsList[[33]]


cov1 <- cov(samples[,c(1,3)])
cov1

cov2 <- cov(samples[,c(2,4)])
cov2

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$removeSamplers(c('a', 'b'))

conf$addSampler(target = c('a[1]', 'b[1]'),
                type = 'RW_block',
                propCov = cov1)


conf$addSampler(target = c('a[2]', 'b[2]'),
                type = 'RW_block',
                propCov = cov2)

conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Rmcmc$samplerFunctions$contentsList[[33]]
Rmcmc$samplerFunctions$contentsList[[33]]$target
Rmcmc$samplerFunctions$contentsList[[33]]$adaptive
Rmcmc$samplerFunctions$contentsList[[33]]$adaptInterval
Rmcmc$samplerFunctions$contentsList[[33]]$propCov
Rmcmc$samplerFunctions$contentsList[[34]]$propCov

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

system.time(samples <- runMCMC(Cmcmc, 10000))

cor(samples)
apply(samples, 2, effectiveSize)
round(samplesSummary(samples), 2)


## compareMCMCs:

library(compareMCMCs)

mi <- list(code = code,
           constants = constants,
           data = data,
           inits = inits)
           
MCMCs <- c('nimble',   ## nimble default
           'nimble_slice', ## slice sampling
           'block1',   ## block (a1,b1), (a2,b2)
           'block2')   ## w/ custom covariance

defs <- list(
    ##block1 = function() {      ## OLD
    block1 = function(model) {   ## correct way
        ## accept the "model" argument
        ## and must return the "conf" object
        conf <- configureMCMC(model)#  "model"
        conf$removeSamplers(c('a', 'b'))
        conf$addSampler(target=c('a[1]','b[1]'),
                        type = 'RW_block')
        conf$addSampler(target=c('a[2]','b[2]'),
                        type = 'RW_block')
        print('in function 1')
        return(conf)
    },
    ##block2 = function() {      ## OLD
    block2 = function(model) {   ## correct way
            conf <- configureMCMC(model)#  "model"
            conf$removeSamplers(c('a', 'b'))
            conf$addSampler(target=c('a[1]','b[1]'),
                            type = 'RW_block',
                            control = list(propCov = cov1))
            conf$addSampler(target=c('a[2]','b[2]'),
                            type = 'RW_block',
                            control = list(propCov = cov2))
            print('yay for CEFE')
            return(conf)
        }
)

ctrl <- list(niter = 20000,
             nburnin = 10000)

comp <- compareMCMCs(modelInfo = mi,
                     MCMCs = MCMCs,
                     nimbleMCMCdefs = defs,
                     MCMCcontrol = ctrl)

make_MCMC_comparison_pages(comp,
                           dir = '~/Desktop/cefe',
                           modelName = 'litter')



## Dipper CJS model:

library(nimble)
load('~/Downloads/dipper_data.Rdata')

code <- nimbleCode({
    ## priors
    phi ~ dunif(0, 1)
    ##p ~ dunif(0, 1)
    for(i in 1:2) {
        p[i] ~ dunif(0, 1)
    }
    one ~ dconstraint(p[1] < p[2])
    pDiff <- p[1] - p[2]
    for(i in 1:N) {
        gender[i] ~ dbern(0.5)
        genderID[i] <- gender[i] + 1
        ## condition on first observations
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(j in (first[i]+1):T) {
            x[i, j] ~ dbern(phi * x[i, j-1])
            ##y[i, j] ~ dbern(  p * x[i, j])
            y[i, j] ~ dbern(p[genderID[i]] * x[i, j])
        }
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]

constants <- list(N = N, T = T, first = first)

data <- list(y = sightings, one = 1)

xInit <- sightings
xInit[xInit == 0] <- 1

inits <- list(phi = 0.5,
              ##p = 0.5,
              p = c(0.5, 0.6),  ## FIXED
              x = xInit,
              gender = rep(0, N))

Rmodel <- nimbleModel(code, constants, data, inits)



Rmodel$calculate()

Rmodel$gender
Rmodel$genderID

Rmodel$initializeInfo()

Rmodel$calculate('x')
Rmodel$calculate('y')
Rmodel$calculate('phi')
Rmodel$calculate('p')
Rmodel$calculate('gender')

Rmodel$calculate()


conf <- configureMCMC(Rmodel)
conf$printSamplers()

conf$printMonitors()
conf$addMonitors('pDiff', 'gender')

conf$printMonitors()

Rmcmc <- buildMCMC(conf)


Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samplesList <- runMCMC(Cmcmc, niter = 10000,
                       nburnin = 5000,
                       nchains = 3,
                       samples = FALSE,
                       summary = TRUE)


samplesList <- runMCMC(Cmcmc, niter = 10000,
                       nburnin = 5000)
samplesList
round(samplesSummary(samplesList), 2)

library(basicMCMCplots)

chainsPlot(samples)
chainsPlot(samples)
samplesSummary



######################################
######################################
######################################
######        Day 2     ##############
######################################
######################################
######################################




##
## Session 3: nimbleEcology R Package
##

## (first) original dipper model:

library(nimble)
load('~/Downloads/dipper_data.Rdata')


code <- nimbleCode({
    ## priors
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        ## condition on first observations
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(j in (first[i]+1):T) {
            x[i, j] ~ dbern(phi * x[i, j-1])
            y[i, j] ~ dbern(  p * x[i, j])
        }
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]

constants <- list(N = N, T = T, first = first)

data <- list(y = sightings)

xInit <- sightings
xInit[xInit == 0] <- 1

inits <- list(phi = 0.5,
              p = 0.5,
              x = xInit)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## -1111.808

conf <- configureMCMC(Rmodel)

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 20000, nburnin = 10000)

samplesList <- runMCMC(Cmcmc, 10000, nburnin = 5000,
                       nchains = 3,
                       samplesAsCodaMCMC = TRUE)

library(coda)
coda::gelman.diag(samplesList)

samplesPlot(samples)


## now: modify dipper model to use dCJS_ss:

library(nimbleEcology)

codeCJS <- nimbleCode({
    ## priors
    ##phi ~ dunif(0, 1)
    ##for(i in 1:6) {
    ##    logit(phi[i]) <- beta0 + beta1 * covariate[i]
    ##}
    p ~ dunif(0, 1)
    for(i in 1:N) {
        logit(phi[i]) <- beta0 + beta1 * ind_covariate[i]
        y[i, first[i]:T] ~ dCJS_ss(probSurvive = phi[i],
                                   probCapture = p,
                                   len = T-first[i]+1)
        #### condition on first observations
        ##x[i, first[i]] <- 1
        ##y[i, first[i]] <- 1
        ##for(j in (first[i]+1):T) {
        ##    x[i, j] ~ dbern(phi * x[i, j-1])
        ##    y[i, j] ~ dbern(  p * x[i, j])
        ##}
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]

constants <- list(N = N, T = T, first = first)

data <- list(y = sightings)

##xInit <- sightings
##xInit[xInit == 0] <- 1

initsCJS <- list(phi = 0.5, p = 0.5)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## -361.8705

conf <- configureMCMC(Rmodel)
conf$printMonitors()
conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 20000, nburnin = 10000)
samplesCJS <- samples

samplesList <- list(latent = samples,
                    dCJS = samplesCJS)

library(basicMCMCplots)

chainsPlot(samplesList)

chainsSummary(samplesList, jitter = 0.2,
              buffer.right = 1, buffer.left = 1)


## finally: compare latent and dCJS models
## using compareMCMCs()
## renameMCMC(comp, newName, oldName)
## make_MCMC_comparison_pages(comp, dir, modelName)

mi1 <- list(code = code, constants = constants,
            data = data, inits = inits)

comp1 <- compareMCMCs(mi1, MCMCs = 'nimble')

miCJS <- list(code = codeCJS, constants = constants,
              data = data, inits = initsCJS)

comp2 <- compareMCMCs(miCJS, MCMCs = 'nimble')

comp1 <- renameMCMC(comp1, 'latent', 'nimble')
comp2 <- renameMCMC(comp2, 'dCJS', 'nimble')

comp <- c(comp1, comp2)


## argument: pageComponents
## list(timing = FALSE, efficiencySummary = FALSE, 
##             efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, 
##             efficiencyDetails = TRUE, posteriorSummary = TRUE)

pgs <- list(paceSummaryAllParams = FALSE,
            posteriorSummary = FALSE)


make_MCMC_comparison_pages(comp,
                           dir = '~/Desktop/cefe',
                           modelName = 'dipper',
                           pageComponents = pgs)


## Occupancy model:

library(nimble)

load('~/Downloads/occupancy_data.Rdata')

length(vege)
dim(wind)
head(obs, 20)

code <- nimbleCode({
    ## priors for alpha and beta coeffs
    for(i in 1:2) {
        alpha[i] ~ dnorm(0, sd = 10000)
        beta[i]  ~ dnorm(0, sd = 10000)
    }
    for(i in 1:S) {
        logit(phi[i]) <- alpha[1] + alpha[2] * vege[i]
        z[i] ~ dbern(phi[i])  ## true pres/absences states
        for(j in 1:T) {
            logit(p[i,j]) <- beta[1] + beta[2] * wind[i,j]
            y[i,j] ~ dbern(p[i,j] * z[i])  ## likelihood
        }
    }
})

S <- dim(wind)[1]
T <- dim(wind)[2]

constants <- list(S = S, T = T, wind = wind, vege = vege)

data <- list(y = obs)

inits <- list(alpha = rep(0, 2), beta = rep(0, 2),
              z = rep(1,S))

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()   ## -317.776

conf <- configureMCMC(Rmodel)

conf$printMonitors()

conf$printSamplers(byType = TRUE)
conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samples <- runMCMC(Cmcmc, 20000, nburnin = 10000)


samplesPlot(samples)

samplesList <- runMCMC(Cmcmc, 20000, nburnin = 10000,
                       nchains = 3,
                       samplesAsCodaMCMC = TRUE)

gelman.diag(samplesList)

round(samplesSummary(samples), 3)


## using dOcc_v(probOcc, probDetect, len) distribution
##                       probDetect is a vector



codeOcc <- nimbleCode({
    ## priors for alpha and beta coeffs
    for(i in 1:2) {
        alpha[i] ~ dnorm(0, sd = 10000)
        beta[i]  ~ dnorm(0, sd = 10000)
    }
    for(i in 1:S) {
        logit(phi[i]) <- alpha[1] + alpha[2] * vege[i]
        logit(p[i,1:T]) <- beta[1] + beta[2] * wind[i,1:T]
        ##
        y[i, 1:T] ~ dOcc_v(probOcc = phi[i],
                           probDetect = p[i, 1:T],
                           len = T)
        ##z[i] ~ dbern(phi[i])  ## true pres/absences states
        ##for(j in 1:T) {
        ##    logit(p[i,j]) <- beta[1] + beta[2] * wind[i,j]
        ##    ##y[i,j] ~ dbern(p[i,j] * z[i])  ## likelihood
        ##}
    }
})

S <- dim(wind)[1]
T <- dim(wind)[2]

constants <- list(S = S, T = T,
                  wind = wind, vege = vege)

data <- list(y = obs)

inits <- list(alpha = rep(0, 2),
              beta = rep(0, 2))
              ##z = rep(1,S))

Rmodel <- nimbleModel(codeOcc, constants, data, inits)
Rmodel$calculate()   ## -317.776

conf <- configureMCMC(Rmodel)

conf$printMonitors()
conf$printSamplers()

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

samplesOcc <- runMCMC(Cmcmc, 20000, nburnin = 10000)

samplesList <- list(latent = samples,
                    dOcc = samplesOcc)

library(basicMCMCplots)

chainsPlot(samplesList)

chainsSummary(samplesList, jitter = 0.2,
              buffer.right = 1, buffer.left = 1)


SL2 <- runMCMC(Cmcmc, 20000, nburnin = 10000, nchains = 4)

chainsSummary(SL2, buffer.right = 0.4, buffer.left = 0.4,
              scale = TRUE)


##
## Session 4: Programming and Advanced Uses of NIMBLE
##


## using nimbleMCMC:

library(nimble)

load('~/Downloads/dipper_data.Rdata')

code <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(j in (first[i]+1):T) {
            x[i, j] ~ dbern(phi * x[i, j-1])
            y[i, j] ~ dbern(p   * x[i, j])
        }
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]

constants <- list(N = N, T = T, first = first)

data <- list(y = sightings)

xInit <- sightings
xInit[xInit == 0] <- 1

inits <- list(phi = 0.5,
              p = 0.5,
              x = xInit)


samples <- nimbleMCMC(code, constants,
                      data, inits,
                      niter = 10000,
                      nburnin = 5000,
                      nchains = 2,
                      monitors = c('p', 'phi', 'x'),
                      samples = FALSE,
                      summary = TRUE,
                      WAIC = TRUE)

out <- samples

str(out)



## estimating the posterior mode for the 
## marginalized occupancy model


load('~/Downloads/occupancy_data.Rdata')


codeOcc <- nimbleCode({
    for(i in 1:2) {
        alpha[i] ~ dnorm(0, sd = 10000)
        beta[i]  ~ dnorm(0, sd = 10000)
    }
    for(i in 1:S) {
        logit(phi[i]) <- alpha[1] + alpha[2] * vege[i]
        logit(p[i,1:T]) <- beta[1] + beta[2] * wind[i,1:T]
        y[i, 1:T] ~ dOcc_v(probOcc = phi[i],
                           probDetect = p[i, 1:T],
                           len = T)
    }
})

S <- dim(wind)[1]
T <- dim(wind)[2]

constants <- list(S = S, T = T,
                  wind = wind,
                  vege = vege)

data <- list(y = obs)

inits <- list(alpha = rep(0, 2),
              beta  = rep(0, 2))

Rmodel <- nimbleModel(codeOcc, constants, data, inits)
Rmodel$calculate()

Rmodel$alpha
Rmodel$beta <- c(-1, -1)
Rmodel$beta
Rmodel$calculate()

Cmodel <- compileNimble(Rmodel)

Cmodel$beta
Cmodel$beta[1]
Cmodel$beta[2]
Cmodel[['beta[1]']] <- 99
Cmodel[['beta[1]']]
Cmodel[['beta[2]']]
Cmodel[['alpha[2]']]

list$element
list[['element']]


f <- function(par) {  ## par is: (a1, a2, b1, b2)
    Cmodel[['alpha[1]']] <- par[1]
    Cmodel[['alpha[2]']] <- par[2]
    Cmodel[['beta[1]']]  <- par[3]
    Cmodel[['beta[2]']]  <- par[4]
    ll <- Cmodel$calculate()
    return(-ll)
}

## inital values to start the optimization at,
## a function to maximize, or minimize

Cmodel$alpha
Cmodel$beta

f( c(0,0,-1,-1) )

system.time(out <- optim(c(0,0,0,0), f))

out$convergence
round(out$par, 2)
      

alpha1 = 0
alpha2 = 3
beta1 = -2
beta2 = -3

nimbleMCMC(codeOcc, constants, data, inits,
           samples = FALSE, summary = TRUE)



## now, the same things, except
## we'll optimize a nimbleFunction

Rmodel <- nimbleModel(codeOcc, constants, data, inits)

Rnf <- nimbleFunction(
    run = function(x = double(1)) {
        y <- x[1:5] + 10
        z <- c(y, x)
        returnType(double(1))
        return(z)
    }
)

Rnf(1)

Cnf <- compileNimble(Rnf)

Rnf(2)
Cnf(2)

Rnf
Cnf

Rnf(1:10)
Cnf(1:10)

Rmodel$getNodeNames(includeData = FALSE,
                    stochOnly = TRUE)

## you can think of this as
## a "nimbleFunction generator"
## this function will *create*
## executeable nimbleFunctions
nfGen <- nimbleFunction(
    setup = function(model) {
        ## runs ONE TIME ONLY
        ## ANY R CODE is ok here
        nodes <- model$getNodeNames(includeData = FALSE,
                                    stochOnly = TRUE)
        message('we are in the setup function')
        message('the nodes were determined to be:')
        print(nodes)
        x <- 10
    },
    run = function(par = double(1)) {
        ## this must be written ONLY using the
        ## NIMBLE DSL language
        ## ("a reduced subset of R")
        ## This function will get compiled
        ##for(i in 1:length(nodes)) {
        ##    model[[nodes[i]]] <<- par[i]
        ##}
        model[['alpha[1]']] <<- par[1]
        model[['alpha[2]']] <<- par[2]
        model[['beta[1]']] <<- par[3]
        model[['beta[2]']] <<- par[4]
        ll <- model$calculate()
        return(-ll)
        returnType(double())
    },
    methods = list(
        addOne = function() {
            x <<- x + 1
        },
        getX = function() {
            returnType(double())
            return(x)
        },
        setX = function(newX = double()) {
            x <<- newX
        },
        sim = function() {
            model$simulate()
        }
    )
)

Rmax <- nfGen(Rmodel)

Rmax$run( c(0,0,-1,-1) )

Rmodel$beta

Cmodel <- compileNimble(Rmodel)
Cmax <- compileNimble(Rmax, project = Rmodel)


Cmax$run( c(0,0,0,0) )

Cmax$getX()
Cmax$addOne()
Cmax$getX()
Cmax$setX(100)
Cmax$sim()

Cmodel$alpha
Cmodel$beta


Cmax$run( c(0,0,-10,-10) )

Rmodel$beta
Cmodel$beta

system.time(out <- optim( c(0,0,0,0), Cmax$run))

system.time(out2 <- optim( c(0,0,0,0), Rmax$run))

out2


## writing MH sampler for the dipper model:

library(nimble)

load('~/Downloads/dipper_data.Rdata')

code <- nimbleCode({
    ##phi ~ dunif(0, 1)
    ##p ~ dunif(0, 1)
    logit.phi ~ dnorm(0, sd = 1000)
    logit.p ~ dnorm(0, sd = 1000)
    phi <- expit(logit.phi)
    p <- expit(logit.p)
    for(i in 1:N) {
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(j in (first[i]+1):T) {
            x[i, j] ~ dbern(phi * x[i, j-1])
            y[i, j] ~ dbern(p   * x[i, j])
        }
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]

constants <- list(N = N, T = T, first = first)

data <- list(y = sightings)

xInit <- sightings
xInit[xInit == 0] <- 1

inits <- list(logit.phi = 0,
              logit.p = 0,
              x = xInit)

Rmodel <- nimbleModel(code, constants, data, inits)

Rmodel$calculate()

## contains = sampler_BASE,
## setup = function(model, mvSaved, target, control)

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$removeSamplers('logit.p', 'logit.phi')
conf$printSamplers()

my_MH <- nimbleFunction(
    name = 'my_MH',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        scale <- control$scale
    },
    run = function() {
        currentValue <- model[[target]]
        proposalValue <- rnorm(1, currentValue, scale)
        currentLP <- model$getLogProb(calcNodes)
        model[[target]] <<- proposalValue
        propLP <- model$calculate(calcNodes)
        logA <- propLP - currentLP
        u <- runif(1, 0, 1)
        if(u < exp(logA)) {
            ## accept
            copy(from=model, to=mvSaved, nodes=calcNodes,
                 row=1, logProb=TRUE)
        } else {
            ## reject
            copy(from=mvSaved, to=model, nodes=calcNodes,
                 row=1, logProb=TRUE)
        }
    },
    methods = list(
        reset = function() {}
    )
)


Rmodel <- nimbleModel(code, constants, data, inits)
conf <- configureMCMC(Rmodel)
conf$removeSamplers('logit.p', 'logit.phi')

ss <- 100

conf$addSampler(target = 'logit.p', type = 'my_MH', control = list(scale = ss))
conf$addSampler(target = 'logit.phi', type = 'my_MH', control = list(scale = ss))

conf$printSamplers()

conf$resetMonitors()
conf$addMonitors(c('phi', 'p'))

Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project = Rmodel,
                       showCompilerOutput = TRUE)

samples5 <- runMCMC(Cmcmc, 10000, nburnin = 5000)

round(samplesSummary(samples), 3)
round(samplesSummary(samples2), 3)
round(samplesSummary(samples3), 3)
round(samplesSummary(samples4), 3)
round(samplesSummary(samples5), 3)


library(basicMCMCplots)
samplesPlot(samples)
samplesPlot(samples2)
samplesPlot(samples3)
samplesPlot(samples4)
samplesPlot(samples5)

SL <- list(s1 = samples,
           ##s2 = samples2,
           s3 = samples3,
           ##s4 = samples4,
           s5 = samples4)

chainsPlot(SL)
           




