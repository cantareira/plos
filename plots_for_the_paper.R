source("data_preparation.r")
source("functions.R")
require(ggplot2)
require(grid)
library(earlywarnings)


## Create the pomp object with the models
c1y.m1 <- create.pomp.3p(cant.plos.1y)
c1y.m0 <- create.pomp.3p.sc(cant.plos.1y)
## Load the fitted models and simulated plausibility intervals
load("fitted_bsmc.RData")
## English locale for time (months, days, etc)
Sys.setlocale(category="LC_TIME", "en_US.UTF-8")
## Data frame with inflow, outflow and volume
Y <- as.data.frame(cant.plos)
cores <- cm.colors(nrow(Y))
## Common themes and scales
t1 <- theme_bw(base_size=20)+
    theme(axis.title.x = element_text(vjust=-1.1),
          axis.title.y = element_text(vjust=2),
          plot.margin=unit(c(1,1,1.5,1),"line"),
          axis.line=element_line(size=1.1),
          legend.position=c(.9,0.2),
          legend.title=element_text(size=0))

t2 <- theme_bw(base_size=20)+
    theme(axis.title.x = element_text(vjust=-1.1),
          axis.title.y = element_text(vjust=2)
          )

escala <- scale_color_gradient(breaks=c(1000,2000,3000),
                               labels=format(min(time(cant.plos))+c(1000,2000,3000),"%Y"),
                               low=cores[1], high="darkblue")

## Time series of volume, rainfall, inflow and outflow
Yz <- cant.plos %>%
        as.data.frame() %>%
            mutate(inflow.m = runmean(inflow, k=30, align="right")/(24*3600),
                   outflow.m = runmean(outflow, k=30, align="right")/(24*3600)) %>%
                select(v.rel, pluv.m, inflow.m, outflow.m) %>%
                    zoo(time(cant.plos))

pdf("volume-rain-flow-ts.pdf", width=12, height=8)
par(mfcol=c(2,1), mar=c(0,7.2,0,6.5), las=1, oma=c(4,4,.5,.5), tcl=-.25, 
    mgp=c(5,1,0), cex.lab=2, cex.axis=1.8, lwd=2.5)
plot(Yz[,"pluv.m"], type="h", col="gray", axes=FALSE, ylab="")
par(las=0)
mtext("Rainfall (mm)", side=4, line=3.8, cex=2)
par(las=1)
axis(4)
par(new=TRUE)
plot(Yz[,"v.rel"], lwd=3, col="blue", axes=FALSE, ylab="Volume (%)")
axis(2, at=seq(0,100, by=20))
par(new=FALSE)
plot(Yz[,"inflow.m"], type="h", col="gray",
     ylab=expression(paste("Water flow (" , m^3*s^-1 ,")", sep="")),
     ylim=c(0,max(Yz[,c("outflow.m","inflow.m")])))
lines(Yz[,"outflow.m"], lwd=3, col="blue")
dev.off()

## Volume x rainfall observed and theoretical
f1a <- ggplot(Y, aes(pluv.m, v.rel, colour=1:nrow(Y))) +
    geom_point()+
        geom_path()+
            xlab("Mean rainfall in previous 30 days (mm)") +
                ylab("Stored water (% operational volume)")+ escala + t1 +
                    xlim(-0.1,22) +
                        theme(legend.position=c(0.92,0.75))

Yt <- read.table("data_funil.dat", col.names=c("time", "v.rel", "pluv.m"))

f1b <- ggplot(Yt, aes(pluv.m, v.rel)) + geom_path()

f1b.inset <- ggplotGrob(f1b + xlab("") + ylab("") +
                            scale_x_continuous(labels=NULL, breaks=NULL) +
                                scale_y_continuous(labels=NULL, breaks=NULL) +
                                    theme_bw() +
                                        theme(plot.margin=unit(c(-0.5,-0.5,-0.5,-.5),"in")))

f1 <- f1a + annotation_custom(grob = f1b.inset, xmin = 13, xmax = 21, ymin = -8, ymax = 25)
f1

## Outflow x volume
f2 <- Y %>%
    filter(outflow<4e6) %>%
        ggplot(aes(v.rel, outflow/(24*3600), colour=1:length(outflow)))+
            geom_point()+
                xlab("Stored water (% operational volume)") +
                    ylab(expression(paste("Water withdraw (",m^3/s,")",sep="")))+
                        escala
print(f2+t1)

## Inflow/rainfall x volume
f3 <- ggplot(Y, aes((v.rel+29.2)/1.292, inflow/(24*3600*(pluv.m+0.1)), colour=1:nrow(Y)))+
        geom_point()+
            xlab("Stored water (% maximum volume)") +
                ylab(expression(paste("Water inflow / rainfall (",m^3/s.mm,")",sep="")))+
                    coord_trans(y="log10", x="log10") +
                        escala +
                            scale_y_continuous(breaks=c(2,5,10,20,40,100,200)) +
                                scale_x_continuous(breaks=c(5, 10,25,50,100))
print(f3 + t1 + theme (legend.position=c(0.15,0.85)))


########################################################
## Difffusion-drift-jump model for conditional variance
########################################################
cant.ddj <- ddjnonparam_ews2(matrix(as.numeric(cant.plos$v.rel)+29.2))
Y1 <- zoo(c(s2=cant.ddj$S2.t), time(cant.plos))
Y2 <- data.frame(vrel=exp(cant.ddj$avec)-29.2, s2=cant.ddj$S2.vec)
f5a <- ggplot(fortify(Y1), aes(x=Index, y=Y1)) +
    geom_line(size=1.2) +
        ylab("") +
            xlab("") +t1
f5b <- ggplot(Y2, aes(x=vrel, y=s2)) +
    geom_line(size=1.2) +
        ylab("Conditional variance") +
            xlab("Stored water (% operational volume)") +t1

multiplot(f5b, f5a, cols=2)

########################################################
## Fitted stochastic model
########################################################
f6 <- ggplot(aes(x=Index, y=Value/10e6), data=fortify(c1y.m1.fit$obs, melt=TRUE)) +
    geom_line(size=2) +
        xlab("") +
            ylab(expression(paste("Stored water (",10^6*m^3,")",sep=""))) +
                geom_line(data = fortify(c1y.m1.fit$summary), aes(x=Index, y=mean/10e6, ymin=lower/10e6, ymax=upper/10e6),
                  size=1.25, colour="blue") +
                    geom_ribbon(data = fortify(c1y.m1.fit$summary), aes(x=Index, y=mean/10e6, ymin=lower/10e6, ymax=upper/10e6),
                                alpha=0.1, fill="blue") +
                        geom_line(data = fortify(c1y.m0.fit$summary), aes(x=Index, y=mean/10e6, ymin=lower/10e6, ymax=upper/10e6),
                                  size=1.25, colour="orange") +
                            geom_ribbon(data = fortify(c1y.m0.fit$summary), aes(x=Index, y=mean/10e6, ymin=lower/10e6,
                                            ymax=upper/10e6), fill="orange",alpha=0.2)
                        
f6 + t1

    

###################################################
## Saving plot images
###################################################
ggsave("volumeXrain.pdf", plot= f1, width=8, height=6)
ggsave("outflowXvolume.pdf", plot= f2 + t1, width=8, height=8)
ggsave("inflow-rainXvolume.pdf", plot= f3 + t1+ theme (legend.position=c(0.15,0.85)), width=8, height=8)
ggsave("sde-fit.pdf", plot= f6 + t1, width=12, height=8)
## To ensure transparency in the publsihed paper - PLOS asks for eps or tiff but eps does not support transparency
ggsave("sde-fit.tiff", plot= f6 + t1, width=7.5, height=5, type="cairo", compression = "lzw")
pdf("conditional-var.pdf", width=12, height=8)
multiplot(f5b, f5a, cols=2)
dev.off()

