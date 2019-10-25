

##
## litters data and example results
##

library(nimble)

littersCode <- nimbleCode({
  for (i in 1:G) {
     for (j in 1:N) {
        # likelihood (data model)
        r[i,j] ~ dbin(p[i,j], n[i,j])
        # latent process (random effects)
        p[i,j] ~ dbeta(a[i], b[i]) 
     }
     # prior for hyperparameters
     a[i] ~ dgamma(1, .001)
     b[i] ~ dgamma(1, .001)
   }
})
##
G <- 2
N <- 16
n <- matrix(c(13, 12, 12, 11, 9, 10, 
              9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 
              7, 5, 10, 7, 6, 10, 10, 10, 7), nrow = 2)
r <- matrix(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 
     12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), 
            nrow = 2)

##save('n', 'r', file = '~/github/nimble/nimble-cefe-2019/docs/data/litters_data.Rdata')

littersConsts <- list(G = G, N = N, n = n)
littersData <- list(r = r)
littersInits <- list( a = c(2, 2), b=c(2, 2) )

littersModel <- nimbleModel(littersCode, 
          data = littersData, constants = littersConsts, inits = littersInits)
cLittersModel <- compileNimble(littersModel)

littersConf <- configureMCMC(littersModel)
littersConf$printMonitors()
littersConf$printSamplers()
littersConf$addMonitors(c('a', 'b', 'p'))
littersMCMC <- buildMCMC(littersConf)

cLittersMCMC <- compileNimble(littersMCMC, project = littersModel)

set.seed(0)
samplesList <- runMCMC(cLittersMCMC,
                       niter = 15000,
                       nburnin = 5000, nchains = 3)

set.seed(0)
samplesCodaList <- runMCMC(cLittersMCMC,
                           niter = 15000,
                           nburnin = 5000, nchains = 3,
                           samplesAsCodaMCMC = TRUE)

library(coda)
ESS <- coda::effectiveSize(samplesList[[1]])

##save('samplesList', 'samplesCodaList', 'ESS', file = '~/github/nimble/nimble-cefe-2019/docs/data/litters_samples.Rdata')



##
## dipper data and example results
##

## creating dipper data
load('~/github/nimble/nimble-cefe-2019/docs/data/dipperData.RData')
sightings <- y
gender <- c(rep(1,124), rep(2,131))
## only keep the 209 individuals that were seen before occasion 6 ( = T-1)
keepInd <- first < 6
sightings <- sightings[keepInd,]
first <- first[keepInd]
gender <- gender[keepInd]
save('sightings', 'first', 'gender', file = '~/github/nimble/nimble-cefe-2019/docs/data/dipper_data.Rdata')

## saving dipper samples
library(nimble)
load('~/github/nimble/nimble-cefe-2019/docs/data/dipper_data.Rdata')
dipperCode <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(t in (first[i]+1):T) {
            x[i, t] ~ dbern(phi * x[i, t-1])
            y[i, t] ~ dbern(p * x[i, t])
        }
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]
dipperConsts <- list(N = N, T = T, first = first)
dipperData <- list(y = sightings)
xInit <- ifelse(!is.na(sightings), 1, 0)
dipperInits <- list(phi = 0.5, p = 0.5, x = xInit)

Rmodel <- nimbleModel(dipperCode, dipperConsts, dipperData, dipperInits)
Rmodel$calculate()   ## now with 209 individuals = -1111.808  ## -1175.578

conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc

set.seed(0)
samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)


## dipper example 2: gender-dependent p, year-dependent phi

dipperCode2 <- nimbleCode({
    for(i in 1:2) {
        p[i] ~ dunif(0, 1)
    }
    for(i in 1:(T-1)) {
        phi[i] ~ dunif(0, 1)
    }
    for(i in 1:N) {
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(t in (first[i]+1):T) {
            x[i, t] ~ dbern(phi[t-1] * x[i, t-1])
            y[i, t] ~ dbern(p[gender[i]] * x[i, t])
        }
    }
})

dipperConsts2 <- list(N = N, T = T, first = first, gender = gender)
dipperData2 <- list(y = sightings)
xInit <- ifelse(!is.na(sightings), 1, 0)
dipperInits2 <- list(phi = rep(0.5,T-1), p = rep(0.5,2), x = xInit)

Rmodel <- nimbleModel(dipperCode2, dipperConsts2, dipperData2, dipperInits2)
Rmodel$calculate()    ## now with 209 individuals = -1111.808  ## -1175.578

conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc

set.seed(0)
samples2 <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)


## dipper example 3: using dCJS_ss

library(nimbleEcology)

dipperCode3 <- nimbleCode({
    phi ~ dunif(0, 1)
    p ~ dunif(0, 1)
    for(i in 1:N) {
        ##x[i, first[i]] <- 1
        ##y[i, first[i]] <- 1
        ##for(t in first[i]:T) {
        ##    x[i, t] ~ dbern(phi * x[i, t-1])
        ##    y[i, t] ~ dbern(p * x[i, t])
        ##}
        y[i, first[i]:T] ~ dCJS_ss(probSurvive = phi,
                                   probCapture = p,
                                   len = T - first[i] + 1)
    }
})

dipperConsts3 <- list(N = N, T = T, first = first)
dipperData3 <- list(y = sightings)
##xInit <- ifelse(!is.na(sightings), 1, 0)
dipperInits3 <- list(phi = 0.5, p = 0.5)

Rmodel <- nimbleModel(dipperCode3, dipperConsts3, dipperData3, dipperInits3)
Rmodel$calculate()   ## now with 209 individuals = -361.8705

conf <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc

set.seed(0)
samples3 <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000)

save('samples', 'samples2', 'samples3', file = '~/github/nimble/nimble-cefe-2019/docs/data/dipper_samples.RData')

library(basicMCMCplots)
chainsSummary(list(orig=samples, dcjs=samples3), jitter = 0.1, buffer.left=1, buffer.right=1)
chainsPlot(list(orig=samples, dcjs=samples3))


## use MCMC comparisons here!!!


mi1 <- list(code = code, constants = constants,
            data = data, inits = inits)

comp1 <- compareMCMCs(modelInfo = mi1,
                      MCMCs = 'nimble',
                      MCMCcontrol = list(niter = 20000, nburnin = 10000))

mi2 <- list(code = code2, constants = constants2,
            data = data2, inits = inits2)

comp2 <- compareMCMCs(modelInfo = mi2,
                      MCMCs = 'nimble',
                      MCMCcontrol = list(niter = 20000, nburnin = 10000))

library(compareMCMCs)

str(comp1)
comp1 <- renameMCMC(comp1, 'orig', 'nimble')
str(comp1)

str(comp2)
comp2 <- renameMCMC(comp2, 'dCJS', 'orig')
str(comp2)


comp_combined <- c(comp1, comp2)

undebug(make_MCMC_comparison_pages)

make_MCMC_comparison_pages(comp_combined,
                           dir = '~/Desktop/cefe',
                           modelName = 'dipper')





##
## occupancy model data and samples
##

## make the occupancy data
set.seed(1)
S <- 100
T <- 3
y <- matrix(NA, nrow = S, ncol = T)
vegHt <- sort(runif(S, -1, 1))
alpha1 <- 0
alpha2 <- 3
psi <- plogis(alpha1 + alpha2 * vegHt)
z <- rbinom(S, 1, psi)
wind <- array(runif(S * T, -1, 1), dim = c(S, T))
beta1 <- -2
beta2 <- -3
p <- plogis(beta1 + beta2 * wind)
for(j in 1:T) y[,j] <- rbinom(S, z, p[,j])
vege <- vegHt
wind <- wind
obs <- y


##save('vege', 'wind', 'obs', file = '~/github/nimble/nimble-cefe-2019/docs/data/occupancy_data.Rdata')


## fitting the occupancy model, and saving results:

library(nimble)

load('~/Downloads/occupancy_data.Rdata')

occCode <- nimbleCode({
    for(i in 1:2) {
        alpha[i] ~ dunif(-100, 100)
        beta[i] ~ dunif(-100, 100)
    }
    for(i in 1:S) {
        logit(psi[i]) <- alpha[1] + alpha[2] * vege[i]
        z[i] ~ dbern(psi[i])
        for(j in 1:T) {
            logit(p[i,j]) <- beta[1] + beta[2] * wind[i,j]
            y[i,j] ~ dbern(z[i] * p[i,j])
        }
    }
})

S <- dim(obs)[1]
T <- dim(obs)[2]
occConsts <- list(S = S, T = T, wind = wind, vege = vege)
occData <- list(y = obs)
occInits <- list(alpha = rep(0, 2), beta = rep(0, 2), z = rep(1, S))

Rmodel <- nimbleModel(occCode, occConsts, occData, occInits)
Rmodel$calculate()   # -298.4521

conf <- configureMCMC(Rmodel)
conf$printSamplers()
conf$printSamplers(byType = TRUE)
conf$printMonitors()

Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc

set.seed(0)
time <- system.time(samples <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000))[3] / 2


## occupancy model using dOcc:

occCode2 <- nimbleCode({
    for(i in 1:2) {
        alpha[i] ~ dunif(-100, 100)
        beta[i] ~ dunif(-100, 100)
    }
    for(i in 1:S) {
        logit(psi[i]) <- alpha[1] + alpha[2] * vege[i]
        ##z[i] ~ dbern(psi[i])
        ##for(j in 1:T) {
        ##    logit(p[i,j]) <- beta[1] + beta[2] * wind[i,j]
        ##    y[i,j] ~ dbern(z[i] * p[i,j])
        ##}
        logit(p[i,1:T]) <- beta[1] + beta[2] * wind[i,1:T]
        y[i,1:T] ~ dOcc_v(probOcc = psi[i], probDetect = p[i,1:T], len = T)
    }
})

S <- dim(obs)[1]
T <- dim(obs)[2]
occConsts2 <- list(S = S, T = T, wind = wind, vege = vege)
occData2 <- list(y = obs)
occInits2 <- list(alpha = rep(0, 2), beta = rep(0, 2))

Rmodel <- nimbleModel(occCode2, occConsts2, occData2, occInits2)
Rmodel$calculate()   # -149.0409
conf <- configureMCMC(Rmodel)
##conf$removeSamplers('beta')
##conf$addSampler('beta', 'RW_block')
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc

set.seed(0)
time2 <- system.time(samples2 <- runMCMC(Cmcmc, niter = 20000, nburnin = 10000))[3] / 2


##save('samples', 'time', 'samples2', 'time2', file = '~/github/nimble/nimble-cefe-2019/docs/data/occupancy_samples.Rdata')




##
## create data for orchid example
##

load('~/github/nimble/nimble-cefe-2019/docs/data/orchidData.RData')
first <- f
save('y', 'first', file = '~/github/nimble/nimble-cefe-2019/docs/data/orchid_data.Rdata')












## writing MH sampler for the dipper model:

library(nimble)
load('~/Downloads/dipper_data.Rdata')

dipperCode <- nimbleCode({
    logit.p ~ dnorm(0, 0.001)
    logit.phi ~ dnorm(0, 0.001)
    p <- expit(logit.p)
    phi <- expit(logit.phi)
    ##phi ~ dunif(0, 1)
    ##p ~ dunif(0, 1)
    for(i in 1:N) {
        x[i, first[i]] <- 1
        y[i, first[i]] <- 1
        for(t in (first[i]+1):T) {
            x[i, t] ~ dbern(phi * x[i, t-1])
            y[i, t] ~ dbern(p * x[i, t])
        }
    }
})

N <- dim(sightings)[1]
T <- dim(sightings)[2]
dipperConsts <- list(N = N, T = T, first = first)
dipperData <- list(y = sightings)
xInit <- ifelse(!is.na(sightings), 1, 0)
dipperInits <- list(logit.phi = 0, logit.p = 0, x = xInit)

samples <- nimbleMCMC(dipperCode, dipperConsts, dipperData, dipperInits,
                      niter = 10000, nburnin = 5000,
                      monitors = c('p', 'phi'))

my_MH <- nimbleFunction(
    name = 'my_MH',
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        calcNodes <- model$getDependencies(target)
        scale <- control$scale
    },
    run = function() {
        initialLP <- model$getLogProb(calcNodes)
        current <- model[[target]]
        proposal <- rnorm(1, current, scale)
        model[[target]] <<- proposal
        proposalLP <- model$calculate(calcNodes)
        lMHR <- proposalLP - initialLP
        if(runif(1,0,1) < exp(lMHR)) {
            ## accept
            copy(from = model, to = mvSaved, nodes = calcNodes, logProb = TRUE, row = 1)
        } else {
            ## reject
            copy(from = mvSaved, to = model, nodes = calcNodes, logProb = TRUE, row = 1)
        }
    },
    methods = list(
        reset = function() {}
    )
)

scale <- 0.05

Rmodel <- nimbleModel(dipperCode, dipperConsts, dipperData, dipperInits)
conf <- configureMCMC(Rmodel, monitors = c('p', 'phi'))
conf$printSamplers()
conf$printSamplers(byType = TRUE)
conf$removeSamplers(c('logit.p', 'logit.phi'))
conf$addSampler(target = 'logit.p', type = 'my_MH', control = list(scale = scale))
conf$addSampler(target = 'logit.phi', type = 'my_MH', control = list(scale = scale))
conf$printSamplers()
conf$printMonitors()
Rmcmc <- buildMCMC(conf)

out <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
Cmcmc <- out$mcmc

samples2 <- runMCMC(Cmcmc, niter = 10000, nburnin = 5000)

samplesSummary(samples2)

chainsPlot(samples2)




## marginalized model: dHMM orchid model:

library(nimble)
library(nimbleEcology)
load('~/Downloads/orchid_data.Rdata')

orchidHMMcode <- nimbleCode({
    ## Survival
    s ~ dunif(0, 1)
    ## Transitions: gamma priors 
    for (i in 1:3){
        a[i] ~ dgamma(1, 1) 
        psiD[i] <- a[i]/sum(a[1:3]) 
        b[i] ~ dgamma(1, 1) 
        psiV[i] <- b[i]/sum(b[1:3]) 
        c[i] ~ dgamma(1, 1) 
        psiF[i] <- c[i]/sum(c[1:3]) 
    }
    ## Define state-transition matrix
    ps[1,1] <- s * psiV[1]
    ps[2,1] <- s * psiV[2]
    ps[3,1] <- s * psiV[3]
    ps[4,1] <- 1-s
    ps[1,2] <- s * psiF[1]
    ps[2,2] <- s * psiF[2]
    ps[3,2] <- s * psiF[3]
    ps[4,2] <- 1-s
    ps[1,3] <- s * psiD[1]
    ps[2,3] <- s * psiD[2]
    ps[3,3] <- s * psiD[3]
    ps[4,3] <- 1-s
    ps[1,4] <- 0
    ps[2,4] <- 0
    ps[3,4] <- 0
    ps[4,4] <- 1
    ## Likelihood
    for (i in 1:N) {
        y[i, first[i]:T] ~ dHMM(init = init[1:4],
                                probObs = po[1:4, 1:3],
                                probTrans = ps[1:4, 1:4],
                                len = length[i])
    }
})

N <- dim(y)[1]
T = dim(y)[2]

length <- T - first + 1
init <- c(1/3, 1/3, 1/3, 0)
po <- array(c(1,0,0,0,0,1,0,0,0,0,1,1), c(4,3))

orchidHMMconsts <- list(N = N, T = T, first = first,
                        length = length, init = init, po = po)

orchidHMMdata <- list(y = y)

orchidHMMinits <- list(a = rep(1, 3),
                       b = rep(1, 3),
                       c = rep(1, 3),
                       s = 0.5)

Rmodel <- nimbleModel(orchidHMMcode, orchidHMMconsts,
                      orchidHMMdata, orchidHMMinits)

Cmodel <- compileNimble(Rmodel)

objective <- function(par, model) {
    nodes <- model$getNodeNames(stochOnly = TRUE, topOnly = TRUE)
    model[[nodes[1]]] <- expit(par[1])
    for(i in 2:10) {
        model[[nodes[i]]] <- exp(par[i])
    }
    lp <- model$calculate()
    return(-lp)
}

## NOTE: need to increase control list maxit here!
out <- optim(par = rep(0,10), fn = objective, model = Cmodel, control = list(maxit = 10000))

str(out)

par <- out$par
par

values <- c(expit(par[1]), exp(par[2:10]))
names(values) <- nodes
values


