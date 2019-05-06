
# ccbio<-read.csv("ccbio2.csv")
# ccenv<-read.csv("ccenv01.csv")
# ccdat<-cbind(ccenv01, ccbio2)
# year<-1982:2013
# ccdat<-cbind(year, ccdat)
# ccdat.long<-melt(ccdat, id=c("year"), na.rm=FALSE, value.name="value")
# 
# ccbio<-read.csv("ccbio.csv")
# ccenv<-read.csv("ccenv.csv")
# ccdat<-cbind(ccenv, ccbio)
# year<-2003:2012
# ccbio<-cbind(year, ccbio)
# ccbio.long<-melt(ccbio, id=c("year"), na.rm=FALSE, value.name="value")
# ccbio.l<-ccbio.long[,c(1,3,2)]
# 
# ccbio<-read.csv("GOMbio2.csv")
# ccenv<-read.csv("GOMenv01.csv")
# ccdat<-cbind(ccenv, ccbio)
# year<-1992:2010
# ccdat<-cbind(year, ccdat)
# ccdat.long<-melt(ccdat, id=c("year"), na.rm=FALSE, value.name="value")
# 
# ccbio<-read.csv("ccbio2.csv")
# ccenv<-read.csv("ccenv01.csv")
# ccdat<-cbind(ccenv, ccbio)
# year<-1964:2013
# ccdat<-cbind(year, ccdat)
# ccdat.long<-melt(ccdat, id=c("year"), na.rm=FALSE, value.name="value")

rm(list = ls())
###############
# Script Info #
###############
# PURPOSE: MARSS DFA approach on ecological indicators
# AUTHOR: Scott Large 2012 (revision 2015)
# REVIEWED BY:
# VERSION: 0.3
#

######################
# CHANGES/ ADDITIONS #
######################
# Added:
#
# Done:
#
############
# PACKAGES #
############
library(MARSS)
library(reshape2)
library(meboot)
library(mgcv)
library(reshape2)
#

set.seed(627)
#cc <- read.csv("data/CCIEA-RPW.csv")
# Subset area... initally for coastwide
#ccALL <- cc[cc$Coastwide == 1,]
#ccALL$year <- as.numeric(ccALL$year)
#
#ccALL$timeseries <- gsub("\\(", "", ccALL$timeseries)
#ccALL$timeseries <- gsub("\\)", "", ccALL$timeseries)
#ccALL$timeseries <- gsub(" ", "_", ccALL$timeseries)
#ccALL$timeseries <- gsub("_-_", "_", ccALL$timeseries)
# Wide df with columns as variables
# dat.full <- ccdat
# #
# #
# ind.name <- names(ccbio)
# driver <- ccenv
# ind<-ccbio.l
# #
# # transpose data, get rid of YEAR and LANDINGS
# df <- ind[!colnames(ind) %in% c("year")]

GOMenv <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/GOMenv.csv")
EBSenv <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/EBSenv.csv")
ccenv <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/ccenv.csv")
NEUSenv <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/NEUSenv.csv")
GOMbio <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/GOMbio.csv")
EBSbio <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/EBSbio.csv")
ccbio <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/ccbio.csv")
NEUSbio <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/NEUSbio.csv")
PIenv <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/PIenv.csv")
PIbio <- read.csv("~/jtwork/IEAnatThresh/IndDFAResults2016/PIbio.csv")

# dataDir <- "~/git/slarge/rCode/data"
# 
# GOMbio <- read.csv(paste0(dataDir, "/GOMbio.csv"))
# GOMenv <- read.csv(paste0(dataDir, "/GOMenv.csv"))

#if NAs in data
NEUSenvfix<-na.roughfix(NEUSenv)


# GOMbio1 <- GOMbio[,-3)] 
# GOMenv1 <- GOMenv[,-c(2, 14, 16,17,18)]
# ccbio1<-ccbio
# ccenv1<-ccenv[,-c(16, 17, 18)]
NEUSbio1<-NEUSbio[,c(2,3,4,5)]
NEUSenv1<-NEUSenvfix[,-11]
EBSbio1<-EBSbio[,c(1,4,5)]
EBSenv1<-EBSenv[,-14]


dat <- t(EBSbio)
### ERROR HERE ####
# I assume the parentheses is a leftover
# driver<-GOMenv1)
driver <- EBSenv1

#dat<-df
N.ts <- dim(dat)[1]
# get length of time series
TT <- dim(dat)[2] 

# Take the z-scores of data to normalize
Sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
Mean <- apply(dat,1,mean,na.rm=TRUE)
dat.z <- (dat-Mean)*(1/Sigma)
rownames(dat.z) <- colnames(dat)
 #save(dri.data, file = "dfaDATA_ccDRI_v0001.Rdata")

##############
# BASE MODEL #
##############
# set new control params
cntl.list = list(minit=200, maxit=1200, allow.degen=FALSE)
# set up forms of R matrices
levels.R <- c("diagonal and equal",
              "diagonal and unequal",
              "unconstrained")
model.data = data.frame()

# NOTE: this will take a long time to run! For the purpose of the dfaGAM we will
# start with N.ts -1 = 4
N.ts = 5
for(R in levels.R) {
  #     for(m in 1:(N.ts-1)) {
  for(m in 1:(N.ts-1)) {    
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(dat.z, 
                 model=dfa.model, 
                 control=cntl.list,
                 silent = T,
                 form="dfa",
                 z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,                                       
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
    cat("Just finished",m,"hidden trend(s) with a",R,"covariance matrix with no covariate \n") # end m loop
  } # end m loop
} # end R loop
# save(file="GOMDFA_ind_v03a.Rdata",list = c("model.data", ls(pattern="^kemz."))) # end R loop
save(file=paste0( "EBSDFA_ind_v01a.Rdata"),list = c("model.data", ls(pattern="^kemz."))) # end R loop

#load("ccDFA_drivers.Rdata")
# load("NEUSDFA_ind_v03a.RDATA")

load(paste0("EBSDFA_ind_v01a.RDATA"))

##############################################################
# Create table of different structures and numbers of trends #
##############################################################
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)

#model.data$Ak.wt = wt/sum(wt)
#
# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-4]
#
best.model <- model.tbl[1,]
fitname = paste("kemz",best.model$m,best.model$R,sep=".")
best.fit = get(fitname)
#a
H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat
#H.inv = promax(coef(best.fit, type="matrix")$Z)$rotmat
#
# rotate factor loadings
Z.rot = coef(best.fit, type="matrix")$Z %*% H.inv   
# rotate trends
year<-1982:2013 #change this accordingly
trends.rot = solve(H.inv) %*% best.fit$states
ts.trends = cbind(year, t(trends.rot))
colnames(ts.trends)<-c("year", "T1", "T2")
#ts.trends = cbind(year, t(trends.rot))
#
save(dat, dat.z, ts.trends, file = paste0( "EBSdfaTrends_v01a.RDATA"))
write.csv(model.tbl, file= paste0("EBSdfaResults_v01a.csv"), row.names = F)
write.csv(ts.trends, file= paste0( "EBSdfaTrendsData_v01a.csv"))

----------------------------------------------------------------------------------------------------------
  #GAM of Indicator DFA and drivers
  
#  driver<- NEUSenv[,-11]
# year<-1964:2013
### ERROR ###
# You probably want to remove year as a response variable
# myresponse <- ts.trends 
myresponse <- ts.trends[,-1] 
mydriver <- cbind(year, driver)
mydriver<-as.matrix(mydriver)
ts.length <- year

### ERROR ###
# No object called ts.trend, and again, you prob don't want year as a response variable
# ind.name<-colnames(ts.trend)
ind.name <- colnames(ts.trends[,-1])

dri.name<-colnames(driver)
# 
name.grid<-expand.grid(ind.name, dri.name)
table.list<-paste(name.grid$Var1,name.grid$Var2, sep=".")

for (d in 1:length(dri.name)){
  for (i in 1:length(ind.name)){
    driver.d<-subset(mydriver, select=c(dri.name)[d])
    driver.yr<-year
    
    indicator.i<-subset(myresponse, select=c(ind.name)[i]) 
    #indicator.NA<-subset(indicator.i, VALUE != "NA")
    #indicator.i<-subset(indicator.NA, YEAR %in% driver.yr)
    indicator.yr<-year
    
    driver.d<-subset(driver.d, year %in% indicator.yr)
    
    driver.indicator.di<-cbind(year, indicator.i,driver.d)
    colnames(driver.indicator.di)<-c("YEAR",ind.name[i],dri.name[d])
    assign(paste(ind.name[i],dri.name[d], sep="."), driver.indicator.di)
    #table.list.i<-paste(ind.name[i],dri.name[d], sep=".")
    #table.list<-c(table.list, table.list.i)
  }
}
# 
# 
# 
# 
# 
# 
# #prev <-progressBar()
# # 
gam.value<-data.frame()

# Since you use i in a for loop below, this might make things more clear
# tl <- "T1.AMO"
for(tl in table.list) {
  # prev<-progressBar(i/length(table.list),prev)
  # Bad karma to name an object after a function... oops
  tbl<-get(tl)         #
  response.name<-colnames(tbl)[2]
  driver.name<-colnames(tbl)[3] 
  myresponse<-tbl[,2]
  mydriver<-tbl[,3]
  ts.length<-tbl[,1]

  
# myresponse<-ts.trends[,3]
# mydriver<-driver[,3]

#
#####################
# Fit the GAM model #
#####################
sp.len <- 200 # Spline length
nb <- 1000  # Number of bootstrap replicates
gam.mat <- matrix(nrow=sp.len, ncol=nb) ## create a matrix to be filled with bootstrapped splines (200 x 1000)
dif1 <- 1 # First derivative
dif2 <- 2 # Second derivative
#
# GAM #

gam1 <- gam(myresponse ~ s(mydriver, bs = "ts"), method = "GCV.Cp", se = T)
# gam1 uses thin-plate splines and prevents overfitting by using cross validation, run 
# ?gam() for more info on model. 
# summary(gam1) will show the summary statistics for the model
# gam.check(gam1)
# plot(gam1) will show a plot of the model

# GLM #
gam2 <- gam(myresponse ~ mydriver, method = "GCV.Cp", se = T)
# gam2 removes the smoothing function s(), so it is basically a linear model. Check out my paper to see the conditions
# on when to pick a GAM over a GLM... although for your data I imagine that GAM will be best.
# gam.check(gam2)

# Pull out relevant model outputs:
summary.gam1 <- as.data.frame(cbind(summary(gam1)$dev.expl,       # deviance explained
                                    summary(gam1)$edf,            # estimated degrees of freedom
                                    summary(gam1)$sp.criterion,   # GCV score
                                    summary(gam1)$s.pv))          # estimated p-value
colnames(summary.gam1)<-c("dev.expl","edf", "GCV", "p_value")

ind.fit <- list(mydriver = seq(from = min(mydriver), to = max(mydriver), length.out = sp.len))
pred <- predict.gam(gam1, newdata = ind.fit, type = "response", se.fit = T)
#
gam.data <- data.frame(pred$fit,
                       ind.fit$mydriver)
colnames(gam.data) <- c("response","driver")

# Bootstrap and subsequent code in a nutshell...
# 1) take a random sample of the data 1000 times
# 2) take GAM of each bootstrap replicate (br)
# 3) Take the 1st and 2nd derivative of each br
# 4) When 95% of the 1st or 2nd derivative replicates are greater or less than zero, we have a significant
#    trend (1st derivative), or threshold (2nd derivative).
#
# Note on maximum entropy bootstrap:
# Bootstrapping is a nonparametric alternative for assessing the 
# uncertainty of a linear trend in the presence
# of autocorrelation. Maximum entropy bootstrap allows
# us to construct a set of replicates of the original series
# to be used for inference while retaining the temporal
# dependence structure of the original series in the resampling process
#

# GAVIN has strong opinions about this, might consider modifying it #
respboot <- meboot(myresponse, reps = nb, trim = 0.1)$ens 
drivboot <- meboot(mydriver, reps = nb, trim = 0.1)$ens
#
gam.mat <- matrix(ncol = nb, nrow = sp.len)
for(i in 1:nb) {
  resp.i <- respboot[,i]
  driv.i <- drivboot[,i]
  gam.i <- gam(resp.i ~ s(driv.i, bs = "ts"), method = "GCV.Cp", se = T)
  newdata.i <- list(driv.i = seq(from = min(driv.i), 
                                 to = max(driv.i),
                                 length.out = sp.len))
  gam.mat[, i] <- predict.gam(gam.i, newdata = newdata.i, type = "response", se.fit = T)$fit
}
#
#######################
# CI of GAM bootstrap #
#######################
ci <- matrix(nrow= 2, ncol= sp.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci)<-c("lower","upper")
for(j in 1:sp.len) {
  IC <- quantile(gam.mat[j,], c(0.025, 0.975))
  ci[,j]<-rbind(IC[1], IC[2])
}

#############################
# 1st Derivative Estimation #
#############################
dif1.line<-diff(gam.data$response, difference=1) # Actual 1st deriv estimate from original smoother
deriv.matrix.1<-matrix(nrow=sp.len-dif1, ncol=nb) ## create a matrix to for the 1st deriv estimates
for(k in 1:nb) {
  derivi<-gam.mat[,k]
  deriv.response<-diff(derivi, difference=dif1)
  driver.len<-length(mydriver)-dif1
  deriv.object <- cbind(deriv.response,driver.len)
  deriv.matrix.1[,k]<-deriv.object[,1]
}


# CI of 1st derivative 
dif1.len <- sp.len-dif1
ci1 <- matrix(nrow= 2, ncol= dif1.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci1) <- c("lower","upper")
#
for(l in 1:dif1.len) {
  IC<-quantile(deriv.matrix.1[l,], c(0.025, 0.975))
  ci1[,l]<-rbind(IC[1], IC[2])
}
#   
#############################
# 2nd Derivative Estimation #
#############################
dif2.line<-diff(gam.data$response, difference=2) # Actual 2nd deriv estimate from original smoother
deriv.matrix.2<-matrix(nrow=sp.len-dif2, ncol=nb) ## create a matrix to for the 2nd deriv estimates
for(m in 1:nb) {
  derivi<-gam.mat[,m]
  deriv.response<-diff(derivi, difference=dif2)
  driver.len<-length(mydriver)-dif2
  deriv.object <- cbind(deriv.response,driver.len)
  deriv.matrix.2[,m]<-deriv.object[,1]
}

# CI of 2nd derivative 
dif2.len<-sp.len-dif2
ci2<-matrix(nrow= 2, ncol= dif2.len) ## create a matrix to be filled with bootstrapped CI
rownames(ci2)<-c("lower","upper")
#
for(n in 1:dif2.len) {
  IC<-quantile(deriv.matrix.2[n,], c(0.025, 0.975))
  ci2[,n]<-rbind(IC[1], IC[2])
}
#
# CI of response
lower <- min(ci[1,  ])
upper <- max(ci[2,  ])
lower1 <- min(ci1[1,  ])
upper1 <- max(ci1[2,  ])
lower2 <- min(ci2[1,  ])
upper2 <- max(ci2[2,  ])
#
CI.table<-matrix(nrow=2, ncol=3, data=c(upper,
                                        lower,
                                        upper1,
                                        lower1,
                                        upper2,
                                        lower2),
                 dimnames=list(c("UPPER", "LOWER"),
                               c("GAM","FIRST","SECOND")))          
CI.mat <- data.frame(CI.table)
#


###################
# Make some plots #
###################
#
changepts.lower1 <- seq(1, length(dif1.line)) [ci1["lower",  ] > 0]
changepts.upper1 <- seq(1, length(dif1.line)) [ci1["upper",  ] < 0]
changepts.lower2 <- seq(1, length(dif2.line)) [ci2["lower",  ] > 0]
changepts.upper2 <- seq(1, length(dif2.line)) [ci2["upper",  ] < 0]
#
driver1 <- gam.data$driver[-1]
driver2 <- gam.data$driver[-c(1,2)]
#
###########
## PLOTS ##
###########
# Will make a plot in your "figures directory"
der.plot <- paste0( tl, "_EBSv01a.png")
# der.plot<-paste0()
png(file = der.plot,
    width = 170, 
    height = 130, 
    units = "mm",
    res = 600)
# 
layout(matrix(c(1:3), 3,1,  byrow = T),
       widths = lcm(8.45), heights = lcm(c(4.95,3.95, 3.95)),
       respect = F)
#
####################################
# Figure A) pressure-response plot #
####################################
par(mfg= c(1,1), 
    mar=c(0,2.15,0,0.25), 
    mgp=c(0,0.25,0), 
    cex= .75, 
    tck=-0.015, 
    family = "sans")

plot(myresponse ~ mydriver,
     ylim = c(CI.mat$GAM[2], CI.mat$GAM[1]),
     type = "n",
     ylab="",
     xlab="",
     axes = F)
#
axis(1, tck = .015, labels= F)
#
gam.ax <- pretty(c(CI.mat$GAM[2], CI.mat$GAM[1]), 3)
axis(2, gam.ax, cex=.75)
mtext(response.name, side= 2,line=1.25, cex= .75)
#
# Shading 95% CIs
polygon(c(gam.data$driver, rev(gam.data$driver)), 
        c(ci[2,  ], rev(ci[1,  ])),
        col = gray(0.95), 
        border = gray(.80))
#
# GAM Line
lines(gam.data$response ~ gam.data$driver, lty =2,  lwd=2)
#
# PLOT DERIVATIVES
points(gam.data$response[ci1["lower",  ] > 0] ~ gam.data$driver[ci1["lower",  ] > 0], 
       col = gray(.25), pch = 16)
points(gam.data$response[ci1["upper",  ] < 0] ~ gam.data$driver[ci1["upper",  ] < 0], 
       col = gray(.25), pch = 16)
points(gam.data$response[ci2["lower",  ] > 0] ~ gam.data$driver[ci2["lower",  ] > 0], 
       col = gray(.35), pch = 16)
points(gam.data$response[ci2["upper",  ] < 0] ~ gam.data$driver[ci2["upper",  ] < 0], 
       col = gray(.35), pch = 16)
#
# Add original data points
points(myresponse~mydriver, pch= 1, cex=.75, lwd=1)
#
# Add legend
legend("topleft",
       legend="a",
       bty = "n",
       box.lwd=0, 
       box.col="white")
#
# Box it up
box()
########################
# Figure B) s'(X) plot #
########################
par(mfg = c(2,1), 
    mar=c(0 ,2.15,0.15,0.25), 
    mgp=c(0,0.25,0), 
    cex= .75, 
    tck=-0.015, 
    family = "sans")
#
plot(dif1.line ~ driver1,
     ylim =  c(CI.mat$FIRST[2], CI.mat$FIRST[1]),
     type= "n",
     xlab= "",
     ylab= "",
     axes = F)
#
axis(1, tck = .015, labels= F)
#
d1.ax <- pretty(c(-.85,.85)*max(abs(CI.mat$FIRST)), 3)
axis(2, d1.ax, cex=.75)
mtext("s'(X)",side= 2,line=1.25, cex= .75)
#
polygon(c(driver1[changepts.upper1], rev(driver1[changepts.upper1])),
        c(ci1["upper",][changepts.upper1], rev(as.matrix(rep(0, length(changepts.upper1)))[,1])),
        col="black", border=NA)
#

polygon(x = c(driver1[changepts.lower1], rev(driver1[changepts.lower1])),
        # c(ci1["lower",][changepts.lower1], rev(as.matrix(rep(0, length(changepts.lower1)))[,1])),
        c(ci1["lower",][changepts.lower1], rev(rep(0, length(changepts.lower1)))),
        col="black", border=NA)
#
polygon(c(driver1, rev(driver1)), 
        c(ci1["upper",], rev(ci1["lower",])),
        col = gray(0.95), border = gray(.80))
#
abline(h=0, lwd=.5)
lines(dif1.line ~ driver1, lwd=1.5)
#
# Arrow stuff should be updated ad hoc
# arrows(x0=gam.data$driver[changepts.lower1][4], y0= min(ci1[1,])*.75,
#        x1=gam.data$driver[changepts.lower1][4], y1= min(ci1[1,])*.15, 
#        col="black", 
#        length=.10,
#        angle=25,
#        lwd=1)
# 
legend("topleft",
       legend="b",
       bty = "n",
       box.lwd=0, 
       box.col="white")
#
box()
########################
# Figure C) s"(X) plot #
########################
par(mfg = c(3,1), 
    mar=c(2.15,2.15,0.15,0.25), 
    mgp=c(0,0.25,0), 
    cex= .75, 
    tck=-0.015, 
    family = "sans")
#
plot(dif2.line ~ driver2,
     ylim =  c(CI.mat$SECOND[1], CI.mat$SECOND[2]),
     type= "n",
     xlab= "",
     ylab= "",
     axes = F)
#
dri.ax <- pretty(range(gam.data$driver), 3)
axis(1, dri.ax, cex = .75)
mtext(driver.name, side = 1, line = 1.25, cex = .75)
#
d2.ax <- pretty(c(-.85,.85) * max(abs(CI.mat$SECOND)), 3)
axis(2, d2.ax, cex=.75)
mtext('s"(X)', side= 2,line=1.25, cex= .75)
#
polygon(c(driver2[changepts.upper2], rev(driver2[changepts.upper2])),
        c(ci2["upper",][changepts.upper2], rev(rep(0, length(changepts.upper2)))),
        col="black", border=NA)
#
### ERROR ###
# This should be the lower changepoints
# polygon(c(driver2[changepts.lower1], rev(driver2[changepts.lower1])),
#         c(ci2[1,][changepts.upper2], rev(as.matrix(rep(0, length(changepts.upper2)))[,1])),
#         col="black", border=NA)
# 
polygon(c(driver2[changepts.lower2], rev(driver2[changepts.lower2])),
        c(ci2["lower",][changepts.lower2], rev(rep(0, length(changepts.lower2)))),
        col="black", border=NA)
#
polygon(c(driver2, rev(driver2)), 
        c(ci2["upper",], rev(ci2["lower",])),
        col = gray(0.95), border = gray(.80))
#
abline(h=0, lwd=.5)
lines(dif2.line ~ driver2, lwd=1.5)
# 
legend("topleft",
       legend="c",
       bty = "n",
       box.lwd = 0, 
       box.col = "white")
#
box()
dev.off()
####
}

###############
#factor loadings plots
#################################
spp = rownames(dat.z)
minZ = 0.2
ylims = c(-1.1*max(abs(Z.rot)), 1.1*max(abs(Z.rot)))
pdf(paste("loadings.",start,".",end,".pdf"))
par(mfrow=c(3,2), mar=c(2,4,1.5,0.5), oma=c(0.4,1,1,1))
for(i in 1:best.model$m) {
  plot(c(1:N.ts)[abs(Z.rot[,i])>minZ], as.vector(Z.rot[abs(Z.rot[,i])>minZ,i]),
       type="h", lwd=2, xlab="", ylab="", xaxt="n", ylim=ylims, xlim=c(0,N.ts+1))
  for(j in 1:N.ts) {
    if(Z.rot[j,i] > minZ) {text(j, -0.05, spp[j], srt=90, adj=1, cex=0.9)}
    if(Z.rot[j,i] < -minZ) {text(j, 0.05, spp[j], srt=90, adj=0, cex=0.9)}
    abline(h=0, lwd=1, col="gray")
  } # end j loop
  mtext(paste("Factor loadings on trend",i,sep=" "),side=3,line=.5)
} # end i loop
dev.off()