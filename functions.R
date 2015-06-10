require(zoo)
require(pomp)
library(grid)
library(earlywarnings)

## Model with effect of stored volume
## Creates a pomp object for a given data series
create.pomp.3p <- function(zoo.obj){
    ## deterministic Skeleton
    skel4b <- function(t,x,params, covars, ...){
        with(as.list(c(x, params, covars)),{
            -outflow + a0*pluv^a1*V^a3
        }
             )
    }
    ## Stochastic Skeleton in C snipet
    skel4b.simC <-
        "double mu;
    mu = V + (-outflow + a0*pow(pluv,a1)*pow(V,a3))*dt;
    if(mu<0.) mu=0.;
    V = rnorm(mu, sigma*V*sqrt(dt));"
    ## time vector
    t1 <- as.numeric(time(zoo.obj))
    t1 <- t1-min(t1)
    ## pomp object with stochastic model
    pomp.obj <- pomp(
        data=data.frame(times=t1,obs=zoo.obj$v.abs),
        times="times",
        t0=0,
        covar=data.frame(
            times=t1,
            pluv=(zoo.obj$pluv.m + 0.05),
            outflow=zoo.obj$outflow),
        tcovar=1,
        statenames="V",
        obsnames="obs",
        paramnames=c("a0","a1","a3","dp","sigma"),
        covarnames=c("pluv", "outflow"),
        measurement.model = obs~norm(mean=V,sd=dp),
        skeleton=skel4b,
        rprocess=euler.sim(step.fun=Csnippet(skel4b.simC), delta.t=0.1),
        parameter.transform=function(params, ...) exp(params),
        parameter.inv.transform=function(params, ...) log(params)
    )
    return(pomp.obj)
}

## Model with effect of rain only
## Creates a pomp object for
create.pomp.3p.sc <- function(zoo.obj){
    ## deterministic Skeleton
    skel4c <- function(t,x,params, covars, ...){
        with(as.list(c(x, params, covars)),{
                 -outflow + a0*pluv^a1
             }
             )
    }
    ## Stochastic Skeleton in C snipet
    skel4c.simC <-
        "double mu;
    mu = V + (-outflow + a0*pow(pluv,a1))*dt;
    if(mu<0.) mu=0.;
    V = rnorm(mu, sigma*V*sqrt(dt));"
    ## time vector
    t1 <- as.numeric(time(zoo.obj))
    t1 <- t1-min(t1)
    ## pomp object with stochastic model
    pomp.obj <- pomp(
        data=data.frame(times=t1,obs=zoo.obj$v.abs),
        times="times",
        t0=0,
        covar=data.frame(
            times=t1,
            pluv=(zoo.obj$pluv.m+0.1),
            outflow=zoo.obj$outflow),
        tcovar=1,
        statenames="V",
        obsnames="obs",
        paramnames=c("a0","a1","dp","sigma"),
        covarnames=c("pluv", "outflow"),
        measurement.model = obs~norm(mean=V,sd=dp),
        skeleton=skel4c,
        rprocess=euler.sim(step.fun=Csnippet(skel4c.simC), delta.t=0.1),
        parameter.transform=function(params, ...) exp(params),
        parameter.inv.transform=function(params, ...) log(params)
        )
    return(pomp.obj)
}

## Fit the 3-parameter model with bsmc2, with fixed parameter gamma (a3)
## using as priors lognormal distribution with parameters taken from another fit
## (for sequential fits as time series are updated)
## if another fit is provided in bsmc2.obj parameters and lognormal priors are taken from it
bsmc2.fit.3p <- function(pomp.obj, params, rpriors, bsmc2.obj, Np=10000, transform=TRUE){
    if(!missing(bsmc2.obj)){
        params <- coef(bsmc2.obj)
        params["V.0"] <- obs(pomp.obj)[1,1]
        post <- bsmc2.obj@post
        if(!bsmc2.obj@transform)
            post <- log(post)
        rpriors <- function(params, ...){
            params["a0"] <- rlnorm(n=1, mean(post["a0",]), sd(post["a0",]))
            params["a1"] <- rlnorm(n=1, mean(post["a1",]), sd(post["a1",]))
            params["sigma"] <- rlnorm(n=1, mean(post["sigma",]), sd(post["sigma",]))
            params
        }
    }
    tmp <- pomp(pomp.obj,
                params=params,
                rprior=rpriors
                )
    bsmc2(tmp, Np=Np, transform=transform)
}

## Fit the 3-parameter model with bsmc2, with free parameter gamma (a3)
## using as priors lognormal distribution with parameters taken from another fit
## (for sequential fits as time series are updated)
## if another fit is provided in bsmc2.obj parameters and lognormal priors are taken from it
bsmc2.fit.3p.a3 <- function(pomp.obj, params, rpriors, bsmc2.obj, Np=10000, transform=TRUE){
    if(!missing(bsmc2.obj)){
        params <- coef(bsmc2.obj)
        params["V.0"] <- obs(pomp.obj)[1,1]
        post <- bsmc2.obj@post
        if(!bsmc2.obj@transform)
            post <- log(post)
        rpriors <- function(params, ...){
            params["a0"] <- rlnorm(n=1, mean(post["a0",]), sd(post["a0",]))
            params["a1"] <- rlnorm(n=1, mean(post["a1",]), sd(post["a1",]))
            params["a3"] <- rlnorm(n=1, mean(post["a3",]), sd(post["a3",]))
            params["sigma"] <- rlnorm(n=1, mean(post["sigma",]), sd(post["sigma",]))
            params
        }
    }
    tmp <- pomp(pomp.obj,
                params=params,
                rprior=rpriors
                )
    bsmc2(tmp, Np=Np, transform=transform)
}

## A function to forecast the fitted stochatic model starting in a given date,
## using a single set of coeficients or a sample from the posterior distribution of coeficients
## p1 is a pomp object and z1 a zoo object with state variables dated, z2 is a zoo object with covars
## all date in z2 should be included in z1
## sim.times = time points at which to simulate, startin with 0 (which corresponds to start)
res.fc <- function(p1, z1, z2, start=min(time(z1)), end=max(time(z1)), deflu,
                   deflu.conv=24*3600, pluv.factor=1,
                   coefs, V.0, nsamp.coef=1, nsim=1000, sim.times, V.max,
                   bounded.vols=TRUE, keep.sims=FALSE, states=TRUE,...){
    coefs <- as.matrix(coefs)
    if(length(time(z1))!=length(time(p1)))
        stop("time lengths differ in pomp and zoo objects")
    if(nsamp.coef > ncol(coefs))
        warning("number of columns in coefs less than nsamp.coef")
    if(!missing(deflu))
       z2$outflow[time(z2)>=start] <- deflu*deflu.conv
    t1 <- zoo(time(p1), time(z1))
    ##t2 <- merge(t1,z2, all=c(FALSE,TRUE))$t1
    t2 <- as.numeric(time(z2)-min(time(p1)))-as.numeric(min(time(z1)))
    p1 <- pomp(
        p1, 
        covar=data.frame(
            times=t2,
            pluv=(z2$pluv.m + 0.1)*pluv.factor,
            outflow=z2$outflow),
        tcovar=1,
        ...
    )
    if(missing(V.0))
        coefs["V.0",] <- as.numeric(z1[time(z1)==start])
    else
        coefs["V.0",] <- V.0
    j <- sample(1:ncol(coefs),nsamp.coef, replace=TRUE)
    t.start <- as.numeric(min(time(p1))+as.Date(start))-as.numeric(min(time(z1)))
    t.end <- as.numeric(min(time(p1))+as.Date(end))-as.numeric(min(time(z1)))
    if(missing(sim.times))
        sim.times <- seq(t.start, t.end)
    else
        sim.times <- t.start + sim.times
    f1 <- function(cfs) {
        pomp::simulate(p1, params=cfs, times=sim.times,
                       nsim=nsim, obs=!states, states=states, as.data.frame=TRUE, t0=t.start)
    }
    sim <- adply(as.matrix(coefs[,j]), 2, f1)
    if(states)
        sim$V2 <- sim$V
    else
        sim$V2 <- sim$obs
    if(bounded.vols){
        sim$V2[is.na(sim$V2)] <- 0
        sim$V2[sim$V2>V.max] <- V.max
    }
    sim.s <-
        sim %>%
            group_by(time) %>%
                summarise(mean=mean(V2, na.rm=TRUE),
                          lower=quantile(V2, 0.025, na.rm=TRUE),
                          upper=quantile(V2, 0.975, na.rm=TRUE),
                          sd= sd(V2, na.rm=TRUE))
    
    sim.s <- zoo(sim.s[,-1], min(time(z1))+sim.times)
    if(keep.sims)
        return(list(sims=sim, obs=z1, summary=sim.s))
    else
        return(list(obs=z1, summary=sim.s))
}

## Forecast for a period ahead
fc.ahead <- function(p1, z1, z2, deflu, pluv.factor=1, ...){
    res.fc(p1=p1, z1=z1,
           z2=z2,
           deflu=deflu,
           pluv.factor=1,
           start=min(time(z2)),
           end=max(time(z2)),
           ...           
           )
}

## Function to plot forecasts generated by function forecast3p
fc.plot <- function(sim, only.obs=FALSE, ci.poly=TRUE, mean.lines=TRUE, ci.lines=FALSE,
                    cpoly=gray.colors(1, alpha=0.3), ...){
    dots <- list(...)
    if(!"ylim" %in% names(dots))
        dots$ylim <- range(c(range(sim$obs),range(sim$summary[,2:3])))
    if(!"col" %in% names(dots))
        dots$col <- "darkblue"
    if(!"lwd" %in% names(dots))
        dots$lwd <- 3
    if(!"ylab" %in% names(dots))
        dots$ylab <- "Stored volume (m3)"
    if(!"xlab" %in% names(dots))
        dots$xlab <-  "Time"
    do.call(plot, c(list(x=sim$obs),dots))
    if(!only.obs){
        if(ci.poly){
            newx <- as.numeric(time(sim$summary))
            y1 <- as.numeric(sim$summary$lower)
            y2 <- as.numeric(sim$summary$upper)
            polygon(c(rev(newx), newx), c(rev(y1),y2), col = cpoly, border = NA)
        }
        if(ci.lines){
            lines(sim$summary$lower, ...)
            lines(sim$summary$upper, ...)
        }
        if(mean.lines)
            if(!"col" %in% names(list(...)))
                dots$col <- "black"
            do.call(lines, c(list(x=sim$summary$mean), dots))
    }
}

## Function to plot lines from forecast generated by function forecast3p
fc.lines <- function(sim, ci.lines=FALSE, mean.lines=TRUE, ci.poly=TRUE,
                     cpoly=gray.colors(1, alpha=0.3), ...){
    if(ci.poly){
        newx <- as.numeric(time(sim$summary))
        y1 <- as.numeric(sim$summary$lower)
        y2 <- as.numeric(sim$summary$upper)
        polygon(c(rev(newx), newx), c(rev(y1),y2), border = NA, col=cpoly)
    }
    if(ci.lines){
        lines(sim$summary$lower, ...)
        lines(sim$summary$upper, ...)
    }
    if(mean.lines)
        lines(sim$summary$mean, ...)
}

## ddjnonparam_ews with logical argument to plot
ddjnonparam_ews2 <- function (timeseries, bandwidth = 0.6, na = 500, logtransform = TRUE, interpolate = FALSE, plot=FALSE) {
    if(class(timeseries)!="matrix") timeseries <- as.matrix(timeseries)
    ## Acessory function not exporte in original package
    Bandi5 <- function (x0, dx, nx, DT, bw, na, avec) 
        {
            SF <- 1/(bw * sqrt(2 * pi))
            x02 <- x0 * x0
            dx2 <- dx * dx
            dx4 <- dx2 * dx2
            dx6 <- dx2 * dx4
            Kmat <- matrix(0, nrow = na, ncol = nx)
            for (i in 1:(nx)) {
                Kmat[, i] <- SF * exp(-0.5 * (x0[i] - avec) * (x0[i] - 
                                                                   avec)/(bw * bw))
            }
            M1.a <- rep(0, na)
            M2.a <- rep(0, na)
            M4.a <- rep(0, na)
            M6M4r <- rep(0, na)
            mean.a <- rep(0, na)
            SS.a <- rep(0, na)
            for (i in 1:na) {
                Ksum <- sum(Kmat[i, ])
                M1.a[i] <- (1/DT) * sum(Kmat[i, ] * dx)/Ksum
                M2.a[i] <- (1/DT) * sum(Kmat[i, ] * dx2)/Ksum
                M4.a[i] <- (1/DT) * sum(Kmat[i, ] * dx4)/Ksum
                M6.c <- (1/DT) * sum(Kmat[i, ] * dx6)/Ksum
                M6M4r[i] <- M6.c/M4.a[i]
                mean.a[i] <- sum(Kmat[i, ] * x0[2:(nx + 1)])/Ksum
                SS.a[i] <- sum(Kmat[i, ] * x02[2:(nx + 1)])/Ksum
            }
            S2.x <- SS.a - (mean.a * mean.a)
            sigma2.Z <- mean(M6M4r)/(5)
            lamda.Z <- M4.a/(3 * sigma2.Z * sigma2.Z)
            sigma2.dx <- M2.a - (lamda.Z * sigma2.Z)
            diff.a <- ifelse(sigma2.dx > 0, sigma2.dx, 0)
            sigma2.dx <- M2.a
            mu.a <- M1.a
            outlist <- list(mu.a, sigma2.dx, diff.a, sigma2.Z, lamda.Z, 
                            S2.x)
            return(outlist)
        }
    timeseries <- ts(timeseries)
    if (dim(timeseries)[2] == 1) {
        Y = timeseries
        timeindex = 1:dim(timeseries)[1]
    }
    else if (dim(timeseries)[2] == 2) {
        Y <- timeseries[, 2]
        timeindex <- timeseries[, 1]
    }
    else {
        warning("not right format of timeseries input")
    }
    if (interpolate) {
        YY <- approx(timeindex, Y, n = length(Y), method = "linear")
        Y <- YY$y
    }
    else {
        Y <- Y
    }
    if (logtransform) {
        Y <- log(Y + 1)
    }
    Xvec1 <- Y
    Tvec1 <- timeindex
    dXvec1 <- diff(Y)
    DT <- Tvec1[2] - Tvec1[1]
    bw <- bandwidth * sd(as.vector(Xvec1))
    alow <- min(Xvec1)
    ahigh <- max(Xvec1)
    na <- na
    avec <- seq(alow, ahigh, length.out = na)
    nx <- length(dXvec1)
    ParEst <- Bandi5(Xvec1, dXvec1, nx, DT, bw, na, avec)
    Drift.vec <- ParEst[[1]]
    TotVar.dx.vec <- ParEst[[2]]
    Diff2.vec <- ParEst[[3]]
    Sigma2Z <- ParEst[[4]]
    LamdaZ.vec <- ParEst[[5]]
    S2.vec <- ParEst[[6]]
    TotVar.i <- approx(x = avec, y = TotVar.dx.vec, xout = Xvec1)
    TotVar.t <- TotVar.i$y
    Diff2.i <- approx(x = avec, y = Diff2.vec, xout = Xvec1)
    Diff2.t <- Diff2.i$y
    Lamda.i <- approx(x = avec, y = LamdaZ.vec, xout = Xvec1)
    Lamda.t <- Lamda.i$y
    S2.i <- approx(x = avec, y = S2.vec, xout = Xvec1)
    S2.t <- S2.i$y
    if(plot){
        dev.new()
        par(mfrow = c(2, 1), mar = c(3, 3, 2, 2), mgp = c(1.5, 0.5, 
                                                      0), oma = c(1, 1, 1, 1))
        plot(Tvec1, Xvec1, type = "l", col = "black", lwd = 2, xlab = "", 
             ylab = "original data")
        grid()
        plot(Tvec1[1:length(Tvec1) - 1], dXvec1, type = "l", col = "black", 
             lwd = 2, xlab = "time", ylab = "first-diff data")
        grid()
        dev.new()
        par(mfrow = c(2, 2), mar = c(3, 3, 2, 2), cex.axis = 1, cex.lab = 1, 
            mgp = c(2, 1, 0), oma = c(1, 1, 2, 1))
        plot(avec, S2.vec, type = "l", lwd = 1, col = "black", xlab = "a", 
             ylab = "conditional variance")
        plot(avec, TotVar.dx.vec, type = "l", lwd = 1, col = "blue", 
             xlab = "a", ylab = "total variance of dx")
        plot(avec, Diff2.vec, type = "l", lwd = 1, col = "green", 
             xlab = "a", ylab = "diffusion")
        plot(avec, LamdaZ.vec, type = "l", lwd = 1, col = "red", 
             xlab = "a", ylab = "jump intensity")
        mtext("DDJ nonparametrics versus a", side = 3, line = 0.1, 
              outer = TRUE)
        dev.new()
        par(mfrow = c(2, 2), mar = c(3, 3, 2, 2), cex.axis = 1, cex.lab = 1, 
            mgp = c(1.5, 0.5, 0), oma = c(1, 1, 2, 1))
        plot(Tvec1, S2.t, type = "l", lwd = 1, col = "black", xlab = "time", 
             ylab = "conditional variance")
        plot(Tvec1, TotVar.t, type = "l", lwd = 1, col = "blue", 
             xlab = "time", ylab = "total variance of dx")
        plot(Tvec1, Diff2.t, type = "l", lwd = 1, col = "green", 
             xlab = "time", ylab = "diffusion")
        plot(Tvec1, Lamda.t, type = "l", lwd = 1, col = "red", xlab = "time", 
             ylab = "jump intensity")
        mtext("DDJ nonparametrics versus time", side = 3, line = 0.1, 
              outer = TRUE)
    }
    nonpar_x <- data.frame(avec, S2.vec, TotVar.dx.vec, Diff2.vec, 
                           LamdaZ.vec)
    nonpar_t <- data.frame(Tvec1, S2.t, TotVar.t, Diff2.t, Lamda.t)
    return(c(nonpar_x, nonpar_t))
}

### Utility functions ###

## Improved ls function (http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session)
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


## Add an alpha value to a colour (http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html#more)
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {


  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
