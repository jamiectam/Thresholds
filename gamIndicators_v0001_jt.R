###############
# Script Info #
###############
# PURPOSE: Generalized Additive Model of indicators
# AUTHOR: Scott Large 2012
# REVIEWED BY:
# vERSION: 0.2
#

######################
# CHANGES/ ADDITIONS #
######################

# Need to add: 
# 1. autocorrelation function

# Done:


############
# PACKAGEs #
############

#library(xtable)
#library(MARSS)
#library(Vennerable)
#library(strucchange)
#library(manipulate)
#library(diptest)
#library(plotrix)
library(mgcv)
library(time)

#source('LMregression.r') Changepoint


###############################################################
## Data should be: csv, txt, or equivalent file
## Columns: time, driver, response
## LTER data: time in years
## SET WORKING DIRECTORY if needed
## REPLACE "mydata.csv" below with your data file
###############################################################

###################
# DATA MANAGEMENT #
###################
if(Sys.info()[['sysname']] == "Windows") {    
  dirstructure <- "~/"
}
if(Sys.info()[['sysname']] == "Darwin") {
  dirstructure <- "/Volumes/WORK FILES/"
}

jtwork <- paste0(dirstructure, "jtwork/IEAnatThresh/data/")
figure.dir <- paste0(dirstructure, "jtwork/IEAnatThresh/data/")
model.dir <- paste0(dirstructure, "jtwork/IEAnatThresh/data/")

setwd(jtwork) 

jtwork<-"/Libraries/Documents/jtwork/IEAnatThresh/data"
#slWork<-"/Users/Scott/NOAA/Work Files/Manuscripts/Thresholds/data"
#figure.hom<-"/Users/Scott/NOAA/Work Files/EcoAP/figures/"
#figure.dir="/Users/slarge/Documents/Manuscripts/Thresholds/figures/"
figure.dir="/Libraries/Documents/jtwork/IEAnatThresh/data"
model.dir="/Libraries/Documents/jtwork/IEAnatThresh/analysis"

setwd(jtwork)


indicator <- read.csv("indicator_test.csv", header=TRUE)
driver <- read.csv("driver_test.csv", header=TRUE)
ts.length<-c(1983:2013)

# Remove Cteno,zoo, limit years 
indicator<-indicator[indicator$YEAR %in% ts.length,]
#                      &
#                      indicator$INDICATOR!= "cteno" &
#                      indicator$INDICATOR!= "zoo",]

# Remove NAO_W, NAO_a, limit years
driver<-driver[driver$YEAR %in% ts.length,]
#                &
#                driver$DRIVER!= "NAO_w" &
#                driver$DRIVER!= "NAO_a",]

ind.name<-unique(as.character(indicator$INDICATOR))
dri.name<-unique(as.character(driver$DRIVER))

name.grid<-expand.grid(ind.name, dri.name)
table.list<-paste(name.grid$Var1,name.grid$Var2, sep=".")

for (d in 1:length(dri.name)){
  for (i in 1:length(ind.name)){
    driver.d<-subset(driver, DRIVER %in% dri.name[d])
    driver.yr<-driver.d$YEAR
    
    indicator.i<-subset(indicator, INDICATOR %in% ind.name[i]) 
    indicator.NA<-subset(indicator.i, VALUE != "NA")
    indicator.i<-subset(indicator.NA, YEAR %in% driver.yr)
    indicator.yr<-indicator.i$YEAR
    
    driver.d<-subset(driver.d, YEAR %in% indicator.yr)
    
    driver.indicator.di<-cbind(driver.d$YEAR,indicator.i$VALUE,driver.d$VALUE)
    colnames(driver.indicator.di)<-c("YEAR",ind.name[i],dri.name[d])
    assign(paste(ind.name[i],dri.name[d], sep="."), driver.indicator.di)
    #table.list.i<-paste(ind.name[i],dri.name[d], sep=".")
    #table.list<-c(table.list, table.list.i)
  }
}






#prev <-progressBar()
# 
gam.value<-data.frame()
for (i in 1:length(table.list)) {
 # prev<-progressBar(i/length(table.list),prev)
  table<-get(table.list[i])         #
  response.name<-colnames(table)[2]
  driver.name<-colnames(table)[3] 
  myresponse<-table[,2]
  mydriver<-table[,3]
  ts.length<-table[,1]
                                }
  #####################
  # Fit the GAM model #
  #####################

  set.seed(1122)
  sp.len<-200 # Spline length
  nb <- 1000  # Bootstrap length
  thresh <- matrix(nrow=sp.len, ncol=nb) ## create a matrix to be filled with bootstrapped splines (200 x 1000)
  dif1<- 1 ## Adjust which difference to take
  dif2<- 2
  
  
  gam1<-gam(myresponse~s(mydriver, bs="ts"), method="GCV.Cp", se=T)
  gam2<-gam(myresponse~mydriver, method="GCV.Cp", se=T)
       
  #gam3<-gam(log(myresponse)~s(mydriver^2, bs="ts"), method="GCV.Cp", se=T)
  #logGam1<-gam(log(myresponse)~s(mydriver, bs="ts"), method="GCV.Cp", se=T)
  
  
  assign(paste(table.list[i], "GAM", sep="."), gam1)
  assign(paste(table.list[i], "LIN", sep="."), gam2)
  #assign(paste(table.list[i], "log", sep="."), logGam1)
  
  
  #value.i<-cbind(table.list[3],summary(gam1)$dev.expl,summary(gam1)$s.table)
  #value.i<-cbind(summary(gam1)$dev.expl, summary(gam1)$s.table)
  value.i<-as.data.frame(cbind(summary(gam1)$dev.expl,summary(gam1)$edf, summary(gam1)$sp.criterion, summary(gam1)$s.pv))   
  value.i<-cbind(table.list[i], "Smoother", as.data.frame(value.i))
  colnames(value.i)<-c("data","model","dev.expl","edf", "GCV", "p_value")
  
  lin.i<-as.data.frame(cbind(summary(gam2)$dev.expl,NA, summary(gam2)$sp.criterion, summary(gam2)$p.pv[2]))
  lin.i<-cbind(table.list[i], "Linear", lin.i)
  colnames(lin.i)<-c("data","model","dev.expl","edf", "GCV", "p_value")

  gam.value<-rbind(gam.value,value.i, lin.i)

    

     # Collect the GAM models that are significant... or at least work)    
na.sig.list<-as.character(gam.value$data[is.na(gam.value$p_value)])
sig.data<-gam.value[!is.na(gam.value$p_value),]
sig.smoother<-sig.data[sig.data$model=="Smoother",]


# Only excludes GAM models that yield <NA> results
# sig.list<-paste(gam.value$data[!is.na(gam.value$p_value)], "GAM", sep=".")
sig.list<-paste(sig.smoother$data[sig.smoother$p_value<0.05], "GAM", sep=".")
# 
# smooth.list<-c("b_scul.LANDINGS.GAM",
#                "b_total.AMO_a.GAM",
#                "b_total.SST.GAM",
#                "b_total.LANDINGS.GAM",
#                "length_mean.LANDINGS.GAM",
#                "length_mean.SST.GAM",
#                "pd_ratio.LANDINGS.GAM",
#                "richness.AMO_a.GAM",  
#                "richness.LANDINGS.GAM",
#                "richness.SST.GAM",
#                "consum.SST.GAM")        

# Check models
for (g in 1:length(sig.list)) {
#  
  resid.plot<-paste(model.dir, sig.list[g],"_rsd",".pdf", sep="")
  pdf(file=resid.plot)  
  plot.gam((get(sig.list[g])), residuals=T, pages=1, pch=19)
  dev.off()
  check.plot<-paste(model.dir, sig.list[g],"_chk",".pdf", sep="")
  pdf(file=check.plot)
  gam.check(get(sig.list[g]), old.style=T)
  dev.off()
                }
  
#write.csv(gam.value, "GAM_results_v001.csv", row.names=F)


# Get the spline for the significant GAM models

for (g in 1:length(sig.list)) {
  gam.table<-get(sig.list[g]) #
  gam.driver<-as.data.frame(get(gsub(".GAM","",sig.list[g])))[,3] #
  pred<-predict.gam(gam.table, se.fit=T)
  #pred <- predict.gam(gam1, se.fit=T)
  fds <- cbind(pred$fit,gam.driver,pred$se.fit)
  fds2 <- order(fds[,2])
  response <- spline(pred$fit[fds2]~gam.driver[fds2],n=sp.len)
  gam.data <- as.data.frame(cbind(response$y,response$x))
  colnames(gam.data)<-c("response","driver")
  assign(paste(sig.list[g], "DATA", sep="."), gam.data)
}

#Bootstrap
#prev <-progressBar()
 for (g in 1:length(sig.list)) {
   #prev<-progressBar(g/length(sig.list),prev)
   sig<-get(paste(sig.list[g], "DATA", sep="."))#
   boot.data<-as.data.frame(get(gsub(".GAM","",sig.list[g])))#
   response.g<-boot.data[,2]
   driver.g<-boot.data[,3]
   year.g<-boot.data[,1]
    for (b in 1:nb) {
      bootsample <- sample(1:length(year.g), replace=T)
      myresponsei <- response.g[bootsample]
      mydriveri <- driver.g[bootsample]
      #######Computing threshold on new sample #######
      gam.object <- gam(myresponsei~s(mydriveri, bs="tp"),method="GCV.Cp", na.action='na.omit')
      prediction <- predict.gam(gam.object, se.fit=T)
      fds <- cbind(prediction$fit,mydriveri,prediction$se.fit,0,0)
      fds2 <- order(fds[,2])
      titu <- spline(prediction$fit[fds2]~mydriveri[fds2],n=sp.len)
      data <- cbind(titu$y,titu$x)
      thresh[,b]<- titu$y
      assign(paste(sig.list[g],"MATRIX",sep="."), thresh)
    }


    #######################
    # CI of GAM bootstrap #
    #######################
    
    ci<-matrix(nrow= 2, ncol= sp.len) ## create a matrix to be filled with bootstrapped CI
    rownames(ci)<-c("lower","upper")
    for (gam.ci in 1:sp.len) {
      IC<-quantile(thresh[gam.ci,], c(0.025, 0.975))
      ci[,gam.ci]<-rbind(IC[1], IC[2])
      assign(paste(sig.list[g],"CI", sep="."), ci)
    }

    #############################
    # 1st Derivative Estimation #
    #############################
    dif1.line<-diff(gam.data$response, difference=1) # Actual 1st deriv estimate from original smoother
    deriv.matrix.1<-matrix(nrow=sp.len-dif1, ncol=nb) ## create a matrix to for the 1st deriv estimates
    for (first in 1:nb) {
      derivi<-thresh[,first]
      deriv.response<-diff(derivi, difference=dif1)
      driver.len<-length(mydriver)-dif1
      deriv.object <- cbind(deriv.response,driver.len)
      deriv.matrix.1[,first]<-deriv.object[,1]
      assign(paste(sig.list[g],"MATRIX.1",sep="."), deriv.matrix.1)
      assign(paste(sig.list[g],"dif1.line", sep="."), dif1.line)
    }
  
   
    # CI of 1st derivative 
    dif1.len<-sp.len-dif1
    ci1<-matrix(nrow= 2, ncol= dif1.len) ## create a matrix to be filled with bootstrapped CI
    rownames(ci1)<-c("lower","upper")
   ##
    #ci1.90<-matrix(nrow= 2, ncol= dif1.len)
    #rownames(ci1.90)<-c("lower","upper")
    for (first.ci in 1:dif1.len) {
      IC<-quantile(deriv.matrix.1[first.ci,], c(0.025, 0.975))
      ci1[,first.ci]<-rbind(IC[1], IC[2])
      assign(paste(sig.list[g], "CI.1", sep="."), ci1)
    # 90% CI of 1st derivative
     # IC.90<-quantile(deriv.matrix.1[first.ci,], c(0.05, 0.95))
      #ci1.90[,first.ci]<-rbind(IC.90[1], IC.90[2])
      #assign(paste(sig.list[g], "CI.1.90", sep="."), ci1.90)
    }
    
    #############################
    # 2nd Derivative Estimation #
    #############################
    dif2.line<-diff(gam.data$response, difference=2) # Actual 2nd deriv estimate from original smoother
    deriv.matrix.2<-matrix(nrow=sp.len-dif2, ncol=nb) ## create a matrix to for the 2nd deriv estimates
    for (second in 1:nb) {
      derivi<-thresh[,second]
      deriv.response<-diff(derivi, difference=dif2)
      driver.len<-length(mydriver)-dif2
      deriv.object <- cbind(deriv.response,driver.len)
      deriv.matrix.2[,second]<-deriv.object[,1]
      assign(paste(sig.list[g], "MATRIX.2",sep="."), deriv.matrix.2)
      assign(paste(sig.list[g],"dif2.line", sep="."), dif2.line)
    }
    
    # CI of 2nd derivative 
    dif2.len<-sp.len-dif2
    ci2<-matrix(nrow= 2, ncol= dif2.len) ## create a matrix to be filled with bootstrapped CI
    rownames(ci2)<-c("lower","upper")
   ##
    #ci2.90<-matrix(nrow= 2, ncol= dif2.len)
    #rownames(ci2.90)<-c("lower","upper") 
   for (second.ci in 1:dif2.len) {
      IC<-quantile(deriv.matrix.2[second.ci,], c(0.025, 0.975))
      ci2[,second.ci]<-rbind(IC[1], IC[2])
      assign(paste(sig.list[g], "CI.2",sep="."), ci2)
    # 90% CI of 1st derivative
     # IC.90<-quantile(deriv.matrix.2[second.ci,], c(0.05, 0.95))
     # ci2.90[,first.ci]<-rbind(IC.90[1], IC.90[2])
     # assign(paste(sig.list[g], "CI.2.90", sep="."), ci2.90)
    }
   
   # CI of response
   lower <- min(ci[1,  ])
   upper <- max(ci[2,  ])
   lower1 <- min(ci1[1,  ])
   upper1 <- max(ci1[2,  ])
   lower2 <- min(ci2[1,  ])
   upper2 <- max(ci2[2,  ])
   
   CI.table<-matrix(nrow=2, ncol=3, data=c(upper,
                                           lower,
                                           upper1,
                                           lower1,
                                           upper2,
                                           lower2),
                    dimnames=list(c("UPPER", "LOWER"),
                                  c("GAM","FIRST","SECOND")))          
   assign(paste(sig.list[g],"CI.matrix",sep="."), CI.table)
   
   # Matrix of 0 values to make polygon for derivative graphs, update with length of significant indicators
   #ze.u1<-as.matrix(rep(0, length(changepts.upper1))) 
   #ze.l1<-as.matrix(rep(0, length(changepts.lower1)))
   #ze.u2<-as.matrix(rep(0, length(changepts.upper2))) 
   #ze.l2<-as.matrix(rep(0, length(changepts.lower2)))
      
}   
  

save.image(file="gamIndicators_v0002_jt.Rdata")

###################
# AUTOCORRELATION #
###################
library(mgcv)
library(nlme)

for (g in 1:length(sig.list)){
  acf.plot<-paste(figure.dir, "acf.", sig.list[g],".pdf", sep="")
  pdf(file=acf.plot)
  
  
  acf.dat<-get(gsub(".GAM","",sig.list[g]))
  acf(residuals.gam(
    gamm(acf.dat[,2]~s(acf.dat[,3]), 
         #method="GCV.Cp",
         correlation=corAR1(0.9, form=~acf.dat[,1]))$gam, type="deviance"),  main=paste("ACF plot for ", sig.list[g], sep=""))
  dev.off()
}

gam1<-gam(myresponse~s(mydriver, bs="ts"), method="GCV.Cp", se=T)


b2<-gamm(acf.dat[,2]~s(acf.dat[,3]), 
     #method="GCV.Cp",
     correlation=corAR1(form=~acf.dat[,1]))

summary(b2$lme)
plot(residuals.gam(b2$gam))

plot(b2$gam, scale=F)

plot(acf.dat[,2],acf.dat[,1])


########   
# PLOT #
########
   
thresholds<-as.data.frame(cbind(RESPONSE=NA,
                  DRIVER=NA,
                  NEG_1l= NA,
                  NEG_1u= NA,
                  POS_1l= NA,
                  POS_1u= NA,
                  NEG_2l= NA,
                  NEG_2u= NA,
                  POS_2l= NA,
                  POS_2u= NA
                  ))


thresholds<-data.frame()
# organize data
for(g in 1:length(sig.list)){
   raw.data<-as.data.frame(get(gsub(".GAM","",sig.list[g])))#
   y.lab.g<-names(raw.data[2])
   x.lab.g<-names(raw.data[3])
   driver.g<-raw.data[,3]
  
   #driver.order<-order(driver.g)
   #driver.g<-driver.g[driver.order]
   response.g<-raw.data[,2]
   #response.g<-response.g[driver.order]
   year.g<-raw.data[,1]
   boot.data<-as.data.frame(get(paste(sig.list[g],"DATA",sep=".")))#
   gam.report<-gam.value[gam.value$data== gsub(".GAM","", sig.list[g]),]#
   
   CI.mat<-as.data.frame(get(paste(sig.list[g],"CI.matrix",sep=".")))
   ci.g<-as.matrix(get(paste(sig.list[g],"CI",sep=".")))
   ci.1.g<-as.matrix(get(paste(sig.list[g],"CI.1",sep=".")))
   ci.2.g<-as.matrix(get(paste(sig.list[g],"CI.2",sep=".")))
   
   dif1.line.g<- diff(get(paste(sig.list[g], ".DATA", sep=""))[,1], difference=1)  
   dif2.line.g<- diff(get(paste(sig.list[g], ".DATA", sep=""))[,1], difference=2)

   changepts.lower1.g <- seq(1, length(dif1.line.g)) [ci.1.g["lower",  ] > 0]
   changepts.upper1.g <- seq(1, length(dif1.line.g)) [ci.1.g["upper",  ] < 0]
   changepts.lower2.g <- seq(1, length(dif2.line.g)) [ci.2.g["lower",  ] > 0]
   changepts.upper2.g <- seq(1, length(dif2.line.g)) [ci.2.g["upper",  ] < 0]

 
   
# Shortened for derivative plots
   dif1.line.dd<- diff(get(paste(sig.list[g], ".DATA", sep=""))[,1], difference=1)
   deriv.matrix.1<-as.matrix(get(gsub(".GAM", ".GAM.MATRIX.1", sig.list[g]),))
   driver1<-boot.data$driver[-1]
   driver2<-boot.data$driver[-c(1,2)]
   
   if (length(changepts.lower1.g)==0){
     POS_1<-c(0,0)
     
   } else {
     POS_1<-range(boot.data$driver[changepts.lower1.g])
   }
         
   if (length(changepts.upper1.g)==0){
     NEG_1<-c(0,0)
     
   } else {
     NEG_1<-range(boot.data$driver[changepts.upper1.g])
   }
   
   if (length(changepts.lower2.g)==0){
     POS_2<-c(0,0)
     
   } else {
     POS_2<-range(boot.data$driver[changepts.lower2.g])
   }
   
   if (length(changepts.upper2.g)==0){
     NEG_2<-c(0,0)
     
   } else {
     NEG_2<-range(boot.data$driver[changepts.upper2.g])
   }   
   
   
   threshx<-cbind(RESPONSE=unlist(strsplit(gsub(".GAM", "", sig.list[g]),"\\."))[1],
                  DRIVER=unlist(strsplit(gsub(".GAM", "", sig.list[g]),"\\."))[2],
                  NEG_1l= as.numeric(NEG_1 [1]),
                  NEG_1u= as.numeric(NEG_1 [2]),
                  POS_1l= as.numeric(POS_1 [1]),
                  POS_1u= as.numeric(POS_1 [2]),
                  NEG_2l= as.numeric(NEG_2 [1]),
                  NEG_2u= as.numeric(NEG_2 [2]),
                  POS_2l= as.numeric(POS_2 [1]),
                  POS_2u= as.numeric(POS_2 [2])
                  )
   thresholds<-rbind(thresholds, threshx)

}
#write.csv(thresholds, "thresholds_v001.csv")
###################################################
# Response-Driver plot w/ highlighted derivatives #
###################################################

   der.plot<-paste(figure.dir, sig.list[g],"_2.pdf", sep="")
#********   pdf(file=der.plot)
   #par(mfrow=c(1,2))
   plot(driver.g,response.g,
        type = "n",
        ylim = c(CI.mat$GAM[2], CI.mat$GAM[1]* 1.05),
        xlab = x.lab.g,
        ylab = y.lab.g,
        bty= "l",
        #main= paste("Inflection points of", y.lab.g, "in response to", x.lab.g, sep=" "))
        main= c("a"))
   
   # Shading 95% CIs
   polygon(c(boot.data$driver, rev(boot.data$driver)), c(ci.g[2,  ], rev(ci.g[1,  ])),
           col = "grey90", border = NA)

   lines(boot.data$response~boot.data$driver, lwd=2)
   
   # Significant 1st derivative points
   points(boot.data$driver[changepts.lower1.g],
          boot.data$response[changepts.lower1.g], 
            lwd=2, pch=19, col="royalblue3")
   points(boot.data$driver[changepts.upper1.g],
          boot.data$response[changepts.upper1.g], 
            lwd=2, pch=19, col="red3")


   #lines(boot.data$driver[changepts.lower1],boot.data$response[changepts.lower1], lwd=5, pch=2, col="red")
   #lines(boot.data$driver[changepts.upper1],boot.data$response[changepts.upper1], lwd=5, pch=2, col="red")
   
   # Significant 2nd derivative inflection points
   points(boot.data$driver[changepts.lower2.g],
          boot.data$response[changepts.lower2.g], 
            lwd=2, pch=19, col="forestgreen")
   points(boot.data$driver[changepts.upper2.g],
          boot.data$response[changepts.upper2.g], 
            lwd=2, pch=19, col="forestgreen")
   
   # Add original data points
   points(response.g~driver.g, pch=20)
   
   # Legend
   #legend("topleft",                       
  #        legend= c("GAM smoother", "Significant 1st Derivative (+)","Significant 1st Derivative (-)",
  #                 "Significant 2nd Derivative",
   #                paste("Deviance explained =", as.numeric(format(gam.report[1,3], digits=4))*100, sep=""),
  #                 paste("p-value =", as.numeric(format(gam.report[1,6], digits=2)), sep="")),  
  #        col= c( "black", "royalblue3", "red3", "forestgreen","", ""),
  #        lty= c(1,1, 1,1,0,0),
  #        lwd= c(2.5,2.5,2.5,2.5,0,0)
          #cex=2.5
        #  )
   
     dev.off()

###################
# PLOT DERIVATIVE #
###################

    #der1.plot<- paste(figure.dir, sig.list[g], "_1stD",".pdf", sep="")
    #pdf(file= der1.plot)
    
   plot(dif1.line.g~driver1,
         type= "n",
         xlab= x.lab.g,
         ylab= "1st derivative of GAM smoother",
         #main= paste("Significant change in 1st derivative \n response of",y.lab.g,"to", x.lab.g, sep=" "),
        main=c("b"),
         ylim= c(CI.mat$FIRST[2], CI.mat$FIRST[1]*1.05))
    
    polygon(c(driver1[changepts.upper1.g], rev(driver1[changepts.upper1.g])),
            c(ci.1.g[2,][changepts.upper1.g], rev(as.matrix(rep(0, length(changepts.upper1.g)))[,1])),
              col="red3", border=NA)

    
    polygon(c(driver1[changepts.lower1.g], rev(driver1[changepts.lower1.g])),
            c(ci.1.g[1,][changepts.lower1.g], rev(as.matrix(rep(0, length(changepts.lower1.g)))[,1])),
              col="royalblue3", border=NA)
              
    polygon(c(driver1, rev(driver1)), 
            c(ci.1.g[2,], rev(ci.1.g[1,])),
              col = "grey90", border = NA)
    abline(h=0)
    lines(dif1.line.g~driver1, lwd=2)
              
              
              
              # add optional arrows to indicate change points
    #change.1 <- driver1[67] # Add year where trend starts
    #change.2 <- driver1[94] # Add year where trend finishes
    # height.y1 <- .15
    # height.y2 <- 0 
    #arrows(change.1, height.y1, change.1, height.y2, col="black", angle=25, lwd=2)
    #arrows(change.2, height.y1, change.2, height.y2, col="black", angle=25, lwd=2)
    dev.off()
}

####################################
# Plot Z-scored response to driver #
####################################

for(g in 1:length(sig.list)){
  z.resp<-as.numeric(scale(
                        get(
                        paste(sig.list[g], "DATA", sep="."))[,1], scale=T, center=T))#
  z.driv<-get(paste(sig.list[g], "DATA", sep="."))[,2]
  z.rel<-as.data.frame(cbind(z.resp[-1],z.driv[-1]))
  colnames(z.rel)<-c("zResp","zDriv")
  z.score<-assign(paste(sig.list[g],"Z.score", sep="."), z.rel)
}
#dif1.line.g<- diff(get(paste(sig.list[g], ".DATA", sep=""))[,1], difference=1)  
#dif2.line.g<- diff(get(paste(sig.list[g], ".DATA", sep=""))[,1], difference=2)

  
zpts.lower1.g<-seq(1, 199) [get(paste(sig.list[15], ".CI.1", sep="")) [1,  ] > 0]
zpts.upper1.g<-seq(1, 199) [get(paste(sig.list[15], ".CI.1", sep="")) [2,  ] < 0]

  
zpts.lower2.g<-seq(1, 198) [get(paste(sig.list[15], ".CI.2", sep="")) [1,  ] > 0]
zpts.upper2.g<-seq(1, 198) [get(paste(sig.list[15], ".CI.2", sep="")) [2,  ] < 0]



# Shortened for derivative plots
land.plot<- paste(figure.dir, "land_plot.pdf",sep="")
pdf(file= land.plot)

plot(pd_ratio.LANDINGS.GAM.Z.score[,1]~
     pd_ratio.LANDINGS.GAM.Z.score[,2],
     type="n",
     ylim=c(-1.25,3),
     ylab="Z-score",
     xlab="LANDINGS (mt)")
 
  lines(b_scul.LANDINGS.GAM.Z.score[,1]~
        b_scul.LANDINGS.GAM.Z.score[,2],
        col="gray30",
        lwd=4,
        lty=1)
  
  lines(pd_ratio.LANDINGS.GAM.Z.score[,1]~
        pd_ratio.LANDINGS.GAM.Z.score[,2],
        col="gray30",
        lwd=4,
        lty=2)

  lines(b_total.LANDINGS.GAM.Z.score[,1]~
    b_total.LANDINGS.GAM.Z.score[,2],
        col="gray60",
        lwd=4,
        lty=1)

  lines(richness.LANDINGS.GAM.Z.score[,1]~
    richness.LANDINGS.GAM.Z.score[,2],
        col="gray60",
        lwd=4,
        lty=2)
rug(x=richness.LANDINGS[,3])

points(b_scul.LANDINGS.GAM.Z.score$zDriv[seq(1, 199)[b_scul.LANDINGS.GAM.CI.1[2,  ] < 0]],
       b_scul.LANDINGS.GAM.Z.score$zResp[seq(1, 199)[b_scul.LANDINGS.GAM.CI.1[2,  ] < 0]], 
       lwd=.5, pch=19, col="red1")

points(pd_ratio.LANDINGS.GAM.Z.score$zDriv[seq(1, 199)[pd_ratio.LANDINGS.GAM.CI.1[2,  ] < 0]],
       pd_ratio.LANDINGS.GAM.Z.score$zResp[seq(1, 199)[pd_ratio.LANDINGS.GAM.CI.1[2,  ] < 0]], 
       lwd=.5, pch=19, col="red1")

points(pd_ratio.LANDINGS.GAM.Z.score$zDriv[seq(1, 199)[pd_ratio.LANDINGS.GAM.CI.2[2,  ] < 0]],
       pd_ratio.LANDINGS.GAM.Z.score$zResp[seq(1, 199)[pd_ratio.LANDINGS.GAM.CI.2[2,  ] < 0]], 
       lwd=.5, pch=19, col="forestgreen")

points(b_total.LANDINGS.GAM.Z.score$zDriv[seq(1, 199)[b_total.LANDINGS.GAM.CI.1[2,  ] < 0]],
       b_total.LANDINGS.GAM.Z.score$zResp[seq(1, 199)[b_total.LANDINGS.GAM.CI.1[2,  ] < 0]], 
       lwd=.5, pch=19, col="red1")

points(b_total.LANDINGS.GAM.Z.score$zDriv[seq(1, 199)[b_total.LANDINGS.GAM.CI.2[1,  ] > 0]],
       b_total.LANDINGS.GAM.Z.score$zResp[seq(1, 199)[b_total.LANDINGS.GAM.CI.2[1,  ] > 0]], 
       lwd=.5, pch=19, col="olivedrab")

points(richness.LANDINGS.GAM.Z.score$zDriv[seq(1, 199)[richness.LANDINGS.GAM.CI.1[2,  ] < 0]],
         richness.LANDINGS.GAM.Z.score$zResp[seq(1, 199)[richness.LANDINGS.GAM.CI.1[2,  ] < 0]], 
         lwd=2, pch=19, col="red1")

legend("topright", 
       legend=c("Sculpin B", "P/D ratio","Total B", "Sp. richness", "s'(X) <0",'s"(X) negative', 's"(X) positive'),
       lty=c(1,2,1,2, 1, 1, 1),
       lwd=c(2, 2, 2, 2, 2, 2, 2),
       col=c("gray30","gray30","gray60","gray60", "red1", "forestgreen", "olivedrab"))

# b_scul gray30 lty=1
#pd_ratio gray 30 lty=2
#b_total gray 60 lty=1
#richness gray 60 lty=2


dev.off()

AMO.plot<-  paste(figure.dir,"AMO_plot.pdf", sep="")
pdf(file= AMO.plot)

plot(b_total.AMO_a.GAM.Z.score[,1]~
  b_total.AMO_a.GAM.Z.score[,2],
     type="n",
     ylim=c(-1.75,1.75),
     ylab="Z-score",
     xlab="AMO")

lines(b_total.AMO_a.GAM.Z.score[,1]~
  b_total.AMO_a.GAM.Z.score[,2],
      col="gray30",
      lwd=4,
      lty=1)

lines(richness.AMO_a.GAM.Z.score[,1]~
  richness.AMO_a.GAM.Z.score[,2],
      col="gray30",
      lwd=4,
      lty=2)

points(b_total.AMO_a.GAM.Z.score$zDriv[seq(1, 199)[b_total.AMO_a.GAM.CI.1[1,  ] > 0]],
       b_total.AMO_a.GAM.Z.score$zResp[seq(1, 199)[b_total.AMO_a.GAM.CI.1[1,  ] > 0]], 
       lwd=.5, pch=19, col="blue1")

points(richness.AMO_a.GAM.Z.score$zDriv[seq(1, 199)[richness.AMO_a.GAM.CI.1[1,  ] > 0]],
       richness.AMO_a.GAM.Z.score$zResp[seq(1, 199)[richness.AMO_a.GAM.CI.1[1,  ] > 0]], 
       lwd=.5, pch=19, col="blue1")

legend("topleft", 
       legend=c("Total B", "Sp. richness", "s'(X) > 0"),
       lty=c(1,2,1),
       lwd=c(2,2,2),
       col=c("gray30","gray30","blue1"))
rug(x=richness.AMO_a[,3])

dev.off()


SST.plot<- paste(figure.dir, "SST_plot.pdf", sep="")
pdf(file= SST.plot)

plot(b_total.SST.GAM.Z.score[,1]~
  b_total.SST.GAM.Z.score[,2],
     type="n",
     ylim=c(-3,3.75),
     ylab="Z-score",
     xlab="SST")

lines(b_total.SST.GAM.Z.score[,1]~
  b_total.SST.GAM.Z.score[,2],
      col="gray30",
      lwd=4,
      lty=1)

lines(length_mean.SST.GAM.Z.score[,1]~
  length_mean.SST.GAM.Z.score[,2],
      col="gray30",
      lwd=4,
      lty=2)

lines(richness.SST.GAM.Z.score[,1]~
  richness.SST.GAM.Z.score[,2],
      col="gray60",
      lwd=4,
      lty=1)

lines(consum.SST.GAM.Z.score[,1]~
  consum.SST.GAM.Z.score[,2],
      col="gray60",
      lwd=4,
      lty=2)

points(b_total.SST.GAM.Z.score$zDriv[seq(1, 199)[b_total.SST.GAM.CI.1[1,  ] > 0]],
       b_total.SST.GAM.Z.score$zResp[seq(1, 199)[b_total.SST.GAM.CI.1[1,  ] > 0]], 
       lwd=.5, pch=19, col="blue1")

points(length_mean.SST.GAM.Z.score$zDriv[seq(1, 199)[length_mean.SST.GAM.CI.1[1,  ] > 0]],
       length_mean.SST.GAM.Z.score$zResp[seq(1, 199)[length_mean.SST.GAM.CI.1[1,  ] > 0]], 
       lwd=.5, pch=19, col="blue1")

points(length_mean.SST.GAM.Z.score$zDriv[seq(1, 199)[length_mean.SST.GAM.CI.2[1,  ] > 0]],
       length_mean.SST.GAM.Z.score$zResp[seq(1, 199)[length_mean.SST.GAM.CI.2[1,  ] > 0]], 
       lwd=.5, pch=19, col="olivedrab")

points(length_mean.SST.GAM.Z.score$zDriv[seq(1, 199)[length_mean.SST.GAM.CI.1[2,  ] < 0]],
       length_mean.SST.GAM.Z.score$zResp[seq(1, 199)[length_mean.SST.GAM.CI.1[2,  ] < 0]], 
       lwd=.5, pch=19, col="red1")

points(richness.SST.GAM.Z.score$zDriv[seq(1, 199)[richness.SST.GAM.CI.1[1,  ] > 0]],
       richness.SST.GAM.Z.score$zResp[seq(1, 199)[richness.SST.GAM.CI.1[1,  ] > 0]], 
       lwd=.5, pch=19, col="blue1")

points(consum.SST.GAM.Z.score$zDriv[seq(1, 199)[consum.SST.GAM.CI.1[1,  ] > 0]],
       consum.SST.GAM.Z.score$zResp[seq(1, 199)[consum.SST.GAM.CI.1[1,  ] > 0]], 
       lwd=.5, pch=19, col="blue1")


rug(x=richness.SST[,3])
legend("top", 
       legend=c("Total B", "Mean Length", "Sp. richness", "Consumption/Landings", "s'(X) > 0", "s'(X) < 0",'s"(X) positive'),
       lty=c(1,2,1,2,1,1,1),
       lwd=c(2,2,2,2,2,2,2),
       col=c("gray30","gray30", "gray60","gray60", "blue1","red1","olivedrab"))


dev.off()


lt.plot<- "landings_time.pdf"
pdf(file= lt.plot)
plot(TL_mean.LANDINGS[,3]~TL_mean.LANDINGS[,1],
     type="l",
     ylab="LANDINGS",
     xlab="YEAR")
abline(h=370000, col="OLIVEDRAB")
abline(h=400000, col="RED")
abline(h=600000, col="RED")
dev.off()

at.plot<- "AMO_time.pdf"
pdf(file= at.plot)
plot(TL_mean.AMO_a[,3]~TL_mean.AMO_a[,1],
     type="l",
     ylab="AMO",
     xlab="YEAR")
abline(h=0.2, col="BLUE")

dev.off()

st.plot<- "SST_time.pdf"
pdf(file= st.plot)

plot(TL_mean.SST[,3]~TL_mean.SST[,1],
     type="l",
     ylab="SST",
     xlab="YEAR")
abline(h=11.5, col="RED")
abline(h=12.75, col="BLUE")

dev.off()

#############################################################################################################################################

#############################################################################################################################################

#gam.ml<-gam(response.g~s(driver.g,bs="tp"),method="REML")
#gam.l<-gam(response.g~driver.g,method="REML")
#gam.gcv<-gam(response.g~s(driver.g, bs="tp"), method="GCV.Cp")
#gam.lg<-gam(response.g~driver.g, method="GCV.Cp")
#gam.p<-gam(response.g~s(driver.g, bs= "ts"), method="REML")
#gam.p.g<-gam(response.g~s(driver.g, bs="ts"), method="GCV.Cp")
#gam.s<-gam(response.g~s(driver.g, bs="tp"), method="GCV.Cp", select=TRUE)

#anova(gam.s, gam.lg)

#plot(gam.s)

#plot(gam.p.g)
#summary(gam.p.g)
#anova(gam.gcv, gam.lg, gam.p.g)

#summary(gam.p)
#anova(gam.p, gam.ml, gam.l)

#summary(gam.lg)
#anova(gam.lg, gam.gcv)

#summary(gam.ml)
#anova(gam.ml,gam.l)
#summary(gam.l)
#plot(gam.ml)
#summary(gam.gcv)
#plot(gam.gcv)
#

    