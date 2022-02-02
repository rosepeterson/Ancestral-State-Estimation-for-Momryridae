require(corHMM)
require(ape)
require(geiger)
require(ggplot2)
require(phytools)
library("MCMCtreeR")
setwd("~/Dropbox/Osteoglossiformes/calibrationtrees_Dec/output/")
#read in pruned time calibrated trees 
ss1ar <- read.tree("ss1abiogeogbears.tree")
ss1cr <- read.tree("ss1cbiogeogbears.tree")
ss2ar <- read.tree("ss2abiogeogbears.tree")
ss2cr <- read.tree("ss2cbiogeogbears.tree")
ss3ar <- read.tree("ss3abiogeogbears.tree")
ss3cr <- read.tree("ss3cbiogeogbears.tree")
ss4ar <- read.tree("ss4abiogeogbears.tree")
ss4cr <- read.tree("ss4cbiogeogbears.tree")
ss5ar <- read.tree("ss5abiogeogbears.tree")
ss5cr <- read.tree("ss5cbiogeogbears.tree")
ss6ar <- read.tree("ss6abiogeogbears.tree")
ss6cr <- read.tree("ss6cbiogeogbears.tree")
ss7ar <- read.tree("ss7abiogeogbears.tree")
ss7cr <- read.tree("ss7cbiogeogbears.tree")
ss8ar <- read.tree("ss8abiogeogbears.tree")
ss8cr <- read.tree("ss8cbiogeogbears.tree")
ss9ar <- read.tree("ss9abiogeogbears.tree")
ss9cr <- read.tree("ss9cbiogeogbears.tree")
ss10ar <- read.tree("ss10abiogeogbears.tree")
ss10ar <- read.tree("ss10cbiogeogbears.tree")
ssacr <- read.tree("ssacbiogeogbears.tree")
trees <- c(ss1ar, ss1cr, ss2ar, ss2cr, ss3ar, ss3cr, ss4ar, ss4cr, ss5ar, ss5cr,
           ss6ar, ss6cr, ss7ar, ss7cr, ss8ar, ss8cr, ss9ar, ss9cr, ss10ar, ss10cr,
           ssacr)

#model testing
fitER <- fitDiscrete(ssacr, x, model = "ER")
fitER
#log-likelihood = -82.770305
#AICc =  167.564563
fitSYM <- fitDiscrete(ss1atree, x, model = "SYM")
fitSYM
#log-likelihood = -70.666150
#AICc = 162.724705
fitARD <- fitDiscrete(ss1atree, x, model = "ARD")
fitARD
#log-likelihood = -66.526858
#AICc = 178.729392
aicc <- c(fitER$opt$aicc, fitSYM$opt$aicc, fitARD$opt$aicc)
aicc

#read in trait data and run simmap on all trees
dat1 <- read.csv("facedata.csv", sep = ",",header=TRUE)#csv of discrete character states for species
facestate = dat1$FACE.STATE..1.normal..2.chin..3.tube..4.barbel..5.barbel.and.tube.
names(facestate) = dat1$Taxa
facestate
facetrees <- make.simmap(trees, facestate, model="SYM", nsim=100)

#summarizes 100 simulation over 21 trees (2100 total trees)
xx<-describe.simmap(facetrees)
xx
layout(matrix(1:210,21,10,byrow=TRUE))
myColors <- c("1" = "#999999", "2" = "#009E73", "3" = "#E69F00", "4" = "Purple", "5" = "#F0E442")

plotSimmap(facetrees[1:10],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[101:110],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[201:210],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[301:310],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[401:410],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[501:510],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[601:610],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[701:710],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[801:810],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[901:910],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1001:1010],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1101:1110],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1201:1210],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1301:1310],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1401:1410],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1501:1510],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1601:1610],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1701:1710],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1801:1810],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[1901:1910],myColors,pts=F,ftype="off", lwd = 0.5)
plotSimmap(facetrees[2001:2010],myColors,pts=F,ftype="off", lwd = 0.5)

### Run ASR on the pruned master Tree using the best AICc model for Figure 4

ssacr <- read.tree("ssacbiogeogbears.tree")
#read in trait data 
traitsrb<-as.vector(read.csv("facedata.csv",row.names=1))
ssarsim <- make.simmap(ssacr, facestate, model="SYM", nsim=100)
PPA <- describe.simmap(ssarsim, message=FALSE)$ace

sodat<-treedata(ssacr, traitsrb, sort=TRUE, warnings=TRUE)
myColors <- c("1" = "#999999", "2" = "#009E73", "3" = "#E69F00", "4" = "Purple", "5" = "#F0E442")
myLabels <- ssarsim$tip.label
myColors

TipColors <- myColors[match(sodat$data[,1], names(myColors), nomatch=1)]#matches trait values to colours using data column 1 in 'dat'
TipColors
names(TipColors)<-ssarsim$tip.label


plot(ssacr, cex = .1)
plotSimmap(ssarsim[1],colors=myColors,lwd=1,ftype="i", fsize = .1,  )
nodelabels(pie=PPA,piecol=statecols.simmap,cex=0.35)
tiplabels(pch=21, col=TipColors, bg = TipColors, cex=0.5)
TipColors <- myColors[match(sodat$data[,1], names(myColors), nomatch=1)]#matches trait values to colours using data column 1 in 'dat'
TipColors
names(TipColors)<-ssarsim$tip.label










XX<-describe.simmap(facetrees)
