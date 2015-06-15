## R codes for paper on catastrophic regime shifts in Sao Paulo Reservoirs
## RM Coutinho, RA Kraenkel adn PI Prado
## Encoding: UTF-8
## prado@ib.usp.br
## Fitting the stochastic models with and without effect of water volume in the reservoir
## The following code refits the model. As the fit is done with probabilistic sampling, the results may differ slightly
## The fit that we used in the paper are in stored in the file fitted_bsmc.RData
## To use this models run the first chunk below (R setup) and then
## c1y.m1 <- create.pomp.3p(cant.plos.1y)
## c1y.m0 <- create.pomp.3p.sc(cant.plos.1y)


###################################################
### R setup
###################################################
source("data_preparation.r")
source("functions.R")
## Maximum capacity of Cantareira
cant.vmax <- 1.2695e9

################################################################################
## 1. Model with effect of rain and volume
################################################################################
## pomp object of the model
c1y.m1 <- create.pomp.3p(cant.plos.1y)
## Starting values: coefs from trajectory matching
guess <- c(a0=9, a1=0.6, a3=0.59, V.0=cant.plos.1y$v.abs[1], dp=3.2e7, sigma=4e-3)
## Checking simulations with starting parameters: bad fit but not nonsense values
fc.plot(res.fc(c1y.m1, cant.plos.1y$v.abs, cant.plos.1y,
               coefs=guess, nsim=100, V.max=cant.vmax))
## BSMC ##
## Priors
c1y.prior <- function(params,...){
    params["a0"] <- runif(n=1, min=0.1, max=20)
    params["a1"] <- runif(n=1, min=0.01,max=0.99)
    params["sigma"] <- runif(n=1, min=0.001,max=0.05)
    params
}
##Fiting with bayesian filter in two steps:
c1y.bsmc <- bsmc2.fit.3p(c1y.m1, params=guess, rpriors=c1y.prior, Np=10000) %>% # 10k thousand iterations with the priors above
    bsmc2.fit.3p(c1y.m1, bsmc2.obj=., Np=30000) ## 30k iterations with lognormal priors from the previous fit
save.image()
## Plausibility intervals of predicted trajectories ##
## With coefficient combinations sampled from posterior distributions
## Trajectories bounded within 0 to maximum volume
c1y.m1.fit <- res.fc(c1y.m1, cant.plos.1y$v.abs, cant.plos.1y,
                           coefs=exp(c1y.bsmc@post),
                     nsamp.coef=5000, nsim=1, V.max=cant.vmax, states=FALSE)
## Quick check
fc.plot(c1y.m1.fit)

################################################################################
## 2. Model with effect of rain only
################################################################################
## pomp object of the model
c1y.m0 <- create.pomp.3p.sc(cant.plos.1y)
## Starting values: coefs from trajectory matching
guess2 <- c(a0=9e5, a1=0.3, V.0=cant.plos.1y$v.abs[1], dp=3.2e7, sigma=4e-3)
## Check simulations with starting parameters
fc.plot(res.fc(c1y.m0, cant.plos.1y$v.abs, cant.plos.1y,
       coefs=guess, nsim=100, V.max=cant.vmax))
## BSMC ##
## Priors
c1y.prior2 <- function(params,...){
    params["a0"] <- runif(n=1, min=1e5, max=1e7)
    params["a1"] <- runif(n=1, min=0.01,max=0.99)
    params["sigma"] <- runif(n=1, min=0.001,max=0.05)
    params
}
##Fiting with bayesian filter in two steps:
c1y.bsmc2 <- bsmc2.fit.3p(c1y.m0, params=guess2, rpriors=c1y.prior2, Np=10000) %>% # 10k iterations with the priors above
    bsmc2.fit.3p(c1y.m0, bsmc2.obj=., Np=30000) # 30k iterations with lognormal priors from the previous fit
save.image()
## Plausibility intervals of predicted trajectories ##
## With coefficient combinations sampled from posterior distributions
## Trajectories bounded within 0 to maximum volume
c1y.m0.fit <- res.fc(c1y.m0, cant.plos.1y$v.abs, cant.plos.1y,
                           coefs=exp(c1y.bsmc2@post),
                     nsamp.coef=5000, nsim=1, V.max=cant.vmax, states=FALSE)
## Quick check
fc.plot(c1y.m0.fit)


## Model comparison: log-likelihoods
LL1 <- c1y.bsmc@log.evidence
LL0 <- c1y.bsmc2@log.evidence
LL1 - LL0

## Saves fitted models and simulated intervals
save(c1y.bsmc, c1y.bsmc2, c1y.m0.fit, c1y.m1.fit, file="fitted_bsmc.Rdata")
