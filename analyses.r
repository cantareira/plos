### R code
## Data set preparation and indicators of critical slow-down
### Encoding: UTF-8

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
## Simulations: worst fit but not nonsense values
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
c1y.bsmc <- bsmc2.fit.3p(c1y.m1, params=guess, rpriors=c1y.prior, Np=10000) %>% ## Ten thousand iterations with the priors above
    bsmc2.fit.3p(c1y.m1, bsmc2.obj=., Np=30000) ## 30k iterations with lognormal priors from the previous fit
save.image()
plot(c1y.bsmc)

## Plausibility intervals of predicted trajectories ##
## With coefficient combinations sampled from posterior distributions
c1y.m1.fit <- res.fc(c1y.m1, cant.plos.1y$v.abs, cant.plos.1y,
                           coefs=exp(c1y.bsmc@post),
                           nsamp.coef=4000, nsim=1, V.max=cant.vmax, states=FALSE)
fc.plot(c1y.m1.fit)

################################################################################
## 2. Model with effect of rain only
################################################################################
## pomp object of the model
c1y.m0 <- create.pomp.3p.sc(cant.plos.1y)
## Starting values: coefs from trajectory matching
guess2 <- c(a0=9e5, a1=0.3, V.0=cant.plos.1y$v.abs[1], dp=3.2e7, sigma=4e-3)
## Simulations: worst fit but not nonsense values
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
c1y.bsmc2 <- bsmc2.fit.3p(c1y.m0, params=guess2, rpriors=c1y.prior2, Np=10000) %>% ## Ten thousand iterations with the priors above
    bsmc2.fit.3p(c1y.m0, bsmc2.obj=., Np=30000) ## 30k iterations with lognormal priors from the previous fit
save.image()
## Plausibility intervals of predicted trajectories ##
## With coefficient combinations sampled from posterior distributions
c1y.m0.fit <- res.fc(c1y.m0, cant.plos.1y$v.abs, cant.plos.1y,
                           coefs=exp(c1y.bsmc2@post),
                           nsamp.coef=4000, nsim=1, V.max=cant.vmax, states=FALSE)
fc.plot(c1y.m0.fit)
## Comparacao dos modelos
aic1 <- -2*c1y.bsmc@log.evidence + 10
aic0 <- -2*c1y.bsmc2@log.evidence + 6
aic0-aic1

################################################################################
## Iterated filter (mif)
################################################################################
save(c1y.bsmc, c1y.bsmc2, file="parallel/needObj.RData")

## Checking results of mif
## Verificando covergencia dos 4 filtros
load("parallel/mif_par3.RData")
mfs <- c(mif.cl.3p[[1]], mif.cl.3p[[2]], mif.cl.3p[[3]], mif.cl.3p[[4]])
plot(mfs)
sapply(mfs, coef)
sapply(mfs, logLik)
