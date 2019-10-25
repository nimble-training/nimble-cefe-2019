
library(nimble)

##
## Session 1: Introduction to NIMBLE
##

littersCode <- nimbleCode({
  for (i in 1:G) {
     for (j in 1:N) {
        # likelihood (data model)
        r[i,j] ~ dbin(p[i,j], n[i,j])
        # latent process (random effects)
        p[i,j] ~ dbeta(a[i], b[i]) 
     }
     # prior for hyperparameters
     # such gamma priors are not generally recommended, but
     # these are the priors from the original example
     a[i] ~ dgamma(1, .001)
     b[i] ~ dgamma(1, .001)
   }
})

n <- matrix(c(13, 12, 12, 11, 9, 10, 
              9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 
              7, 5, 10, 7, 6, 10, 10, 10, 7), nrow = 2)
	      
r <- matrix(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 
     12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), 
     nrow = 2)


## code, constants, data, inits

## NIMBLE model object: littersModel

## compiled model object: cLittersModel





cLittersModel$p
cLittersModel$calculate('a')   # log-prior density
cLittersModel$getLogProb('a')

cLittersModel$a <- c(3, 3)
cLittersModel$getLogProb('a')
cLittersModel$calculate('a')   # log-prior density

set.seed(1)  # so the calculations are reproducible
littersModel$simulate('p')  # simulate from prior
littersModel$p
littersModel$getLogProb('p')  # log prob not yet updated!
littersModel$calculate('p')   # update it
littersModel$getLogProb('p')  # now we're good



##
## Session 2: Improving MCMC Performance, and Dipper Model
##

## slice sampling for litters model:

## block sampling for litters model:

## block sampling with initial covariance matrices:

## compareMCMCs:

library(compareMCMCs)

comp <- compareMCMCs(modelInfo,
                     MCMCs,
                     nimbleMCMCdefs,
                     monitors,
                     MCMCcontrol)

make_MCMC_comparison_pages(comp,
                           dir,
                           modelName)





## Dipper CJS model:

library(nimble)
load('~/Downloads/dipper_data.Rdata')










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


## now: modify dipper model to use dCJS_ss:



## finally: compare latent and dCJS models
## using compareMCMCs()
## renameMCMC(comp, newName, oldName)
## make_MCMC_comparison_pages(comp, dir, modelName)

## argument: pageComponents
## list(timing = FALSE, efficiencySummary = FALSE, 
##             efficiencySummaryAllParams = TRUE, paceSummaryAllParams = TRUE, 
##             efficiencyDetails = TRUE, posteriorSummary = TRUE)




## Occupancy model:

library(nimble)

load('~/Downloads/occupancy_data.Rdata')



## using dOcc_v(probOcc, probDetect, len) distribution






##
## Session 4: Programming and Advanced Uses of NIMBLE
##


## ML estimations for the
## marginalized occupancy model


load('~/Downloads/occupancy_data.Rdata')




## optim(par, function, control = list(fnscale = -1))















