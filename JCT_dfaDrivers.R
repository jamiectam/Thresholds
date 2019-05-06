
# library(xtable)
library(MARSS)
# library(mgcv)
# library(time)
library(reshape2)
# library(car)

##################
# Set up folders #
##################
slWork<-"/Users/slarge/Documents/Manuscripts/MV-Indicator/data"
figure.dir <-"/Users/slarge/Documents/Manuscripts/MV-Indicator/figures/"
source.dir <-"/Users/slarge/Documents/Manuscripts/MV-Indicator/functions/"
# slHome<-"/Users/Scott/Desktop/MV-Indicator/data"
# figure.dir <-"//Users/Scott/Desktop/MV-Indicator/figures/"
setwd(slWork)

# load("gamGradientToyDat.RDATA")
driver <- read.csv("driver_v0004.csv")
indicator <- read.csv("indicator_v0003.csv")

#for R studio
driver<-driver_v0004
indicator<-indicator_v0003


ZOO <- indicator[indicator$INDICATOR == "zoo",]
colnames(ZOO) <- c("YEAR", "VALUE", "DRIVER")
ZOO$DRIVER <- "ZOO"


indicator <- indicator[indicator$INDICATOR != "consum" &
                         indicator$INDICATOR != "zoo",]

driver <- driver[driver$DRIVER != "PRECIP",]
driver <- rbind(driver, ZOO)

dri.wide <- dcast(driver, YEAR ~ DRIVER, value.var = "VALUE")
cc_bio_wide<- dcast(CC_biologic, YEAR~METRIC, value.var="VALUE")
##############

# tf <- dri.wide[complete.cases(dri.wide),]
#tf <- dri.wide[dri.wide$YEAR >= 1964,]
tf <- driver[driver$YEAR >= 1981,]

year <- tf$YEAR
# transpose data, get rid of YEAR and LANDINGS
#df <- tf[,c(-1, -4)]
df <- tf[,c(-1)]

#turn data.frame into matrix                         
dat <- t(df)

N.ts <- dim(dat)[1]
# get length of time series
TT <- dim(dat)[2] 

# Take the z-scores of data to normalize
Sigma <- sqrt(apply(dat,1,var,na.rm=TRUE))
Mean <- apply(dat,1,mean,na.rm=TRUE)
dat.z <- (dat-Mean)*(1/Sigma)
rownames(dat.z) <- colnames(df)

# save(dri.data, file = "dfaDATA_DRI_v0001.Rdata")

########################
# Explore the ENV data #
########################
# dris <- ts(dri.wide[,c(2,4:9)], start= 1950, end = 2011, frequency= 1)
# # plot(dris)
# 
# big <- t(dat.z.scored[,!apply(is.na(dat.z.scored),2, any)])
# big.i <- t(dat[,!apply(is.na(dat),2, any)])
# Using Zuur et al 2010 Methods in Ecology and Evolution protocol
# source("HighstatLib.R") 
# corvif(big.i[,c(11:18)])

################################
# Find optimum model structure #
################################
# if("DFA_COV_model.DRI.LTS.Rdata" %in% dir()) {
#   load("DFA_COV_model.DRI.LTS.Rdata")
#   cov.res = TRUE
#   } else {cov.res = FALSE}
# if("DFA_mCOV_model.DRI.LTS.Rdata" %in% dir()) {
#   load("DFA_mCOV_model.DRI.LTS.Rdata")
#   mcov.res = TRUE
# } else {mcov.res = FALSE}
# if("DFA_model.DRI.LTS.Rdata" %in% dir()) {
#   load("DFA_model.DRI.LTS.Rdata")
#   saved.res = TRUE
# } else {saved.res = FALSE}
# 

##############
# BASE MODEL #
##############
# set new control params
cntl.list = list(minit=200, maxit=1200, allow.degen=FALSE)
# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "unconstrained")
model.data = data.frame()
# r.fit <- data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  #     for(m in 1:(N.ts-1)) {
  for(m in 1:(N.ts-1)) {    
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(dat.z, 
                 model=dfa.model, 
                 control=cntl.list,
                 #control=cntl.list,
                 silent = T,
                 form="dfa",
                 z.score=TRUE)
                 #MCInit=TRUE) #Can't get it to work with MCInit
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,                                       
                                  stringsAsFactors=FALSE))
    r.fit <- rbind(#removed term r.fit here
      data.frame(mean(diag(coef(kemz, type = "matrix")$R)),
                 rbind(diag(coef(kemz, type = "matrix")$R))))
    assign(paste("kemz", m, R, sep="."), kemz)
    cat("Just finished",m,"hidden trend(s) with a",R,"covariance matrix with no covariate \n") # end m loop
  } # end m loop
} # end R loop
save(file="DFA_driversCC.Rdata",list = c("model.data", ls(pattern="^kemz."))) # end R loop
# load("DFA_drivers.Rdata")




##########vignette version

## # set new control params
cntl.list = list(minit=200, maxit=2500, allow.degen=FALSE)
## # set up forms of R matrices
 levels.R = c("diagonal and equal",
              "diagonal and unequal",
              "unconstrained")
model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
 for(R in levels.R) {
     for(m in 1:(N.ts-1)) {
         dfa.model = list(A="zero", R=R, m=m)
         kemz = MARSS(dat.z, model=dfa.model, control=cntl.list, 
             form="dfa", silent=TRUE, z.score=TRUE)
         model.data = rbind(model.data,
                            data.frame(R=R,
                                       m=m,
                                       logLik=kemz$logLik,
                                       K=kemz$num.params,
                                       AICc=kemz$AICc,
                                       stringsAsFactors=FALSE))
         assign(paste("kemz", m, R, sep="."), kemz)
         } # end m loop
     } # end R loop

## # set new control params
cntl.list = list(minit=200, maxit=1200, allow.degen=FALSE)
## # set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "unconstrained")
model.data = data.frame()
# fit lots of models & store results
# NOTE: this will take a long time to run!
for(R in levels.R) {
  for(m in 1:(N.ts-1)) {
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(dat.z, model=dfa.model, control=cntl.list, 
                 form="dfa", silent=TRUE, z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

##############################################################
# Create table of different structures and numbers of trends #
##############################################################
# calculate delta-AICc
model.data$delta.AICc = model.data$AICc - min(model.data$AICc)
# calculate Akaike weights
wt = exp(-0.5*model.data$delta.AICc)
model.data$Ak.wt = wt/sum(wt)

# sort results
model.tbl = model.data[order(model.data$AICc),-4]
# drop AICc from table
# calculate cumulative wts
model.tbl$Ak.wt.cum = cumsum(model.tbl$Ak.wt)
model.tbl = model.tbl[,-4]


best.model <- model.tbl[1,]
fitname = paste("kemz",best.model$m,best.model$R,sep=".")
best.fit = get(fitname)

H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat

#mat1 <- coef(best.fit, type="matrix")$R
#mat2 <- coef(best.fit, type="matrix")$R 
# 
#H.inv = diag(mat2)- diag(mat1)

# rotate factor loadings
Z.rot = coef(best.fit, type="matrix")$Z %*% H.inv   
# rotate trends
trends.rot = solve(H.inv) %*% best.fit$states
ts.trends = cbind(year, t(trends.rot))
colnames(ts.trends) <- c("YEAR", "T1", "T2", "T3")
# colnames(ts.trends) <- c("YEAR", "T1", "T2")

save(driver, indicator, dat, dat.z, ts.trends, file = "dfaTrendsCC_v0001.RDATA")
write.csv(model.tbl, paste0("dfaResults_v0005CC.csv"), row.names = F)


#############
# Load data #
#############
load("dfaTrendsCC_v0001.RDATA")
# land <- driver[driver$DRIVER == "LANDINGS",]
#land <- read.csv("landingsTotal_v006.csv")

#for R studio
land<-landings_cc
indicator<-cc_bio_L

dat.mer <- merge(land[,1:2], ts.trends, by= "YEAR")
colnames(dat.mer) <- c("YEAR", "LANDINGS", "T1", "T2", "T3")

driver <- melt(dat.mer, id.vars= "YEAR")
colnames(driver) <- c("YEAR", "DRIVER", "VALUE")

# Change Landings and total biomass to '000 Tons
driver$VALUE[driver$DRIVER == "LANDINGS"] <- driver$VALUE[driver$DRIVER == "LANDINGS"]/1000
indicator$VALUE[indicator$INDICATOR == "biomass"] <- indicator$VALUE[indicator$INDICATOR == "biomass"]/1000

indicator <- indicator[indicator$INDICATOR != "zoo",]
# unique(indicator$INDICATOR)
# [1] b_total     b_scul      cteno       richness    TL_mean     pd_ratio    plank       length_mean


#################
# Data Analysis #
#################

# Full model
ind.list <- unique(indicator$INDICATOR)
# 
# FM1 <- formula(get(ei.name) ~ s(LANDINGS, k = fm.k3))        # LANDINGS
# FM2 <- formula(get(ei.name) ~ s(T1, k = fm.k3))              # T1
# FM3 <- formula(get(ei.name) ~ s(T2, k = fm.k3))              # T2
# # FM4 <- formula(get(ei.name) ~ s(T1, k = fm.k3) + s(T2, k = fm.k3))        # T1 + T2
# FM4 <- formula(get(ei.name) ~ s(LANDINGS, k = fm.k3) + s(T1, k = fm.k3))  # LANDINGS + T1
# FM5 <- formula(get(ei.name) ~ s(LANDINGS, k = fm.k3) + s(T2, k = fm.k3))  # LANDINGS + T2
# # FM7 <- formula(get(ei.name) ~ s(T1, k = fm.k3) + s(T2, k = fm.k3))        # LANDINGS + T1
# FM.list <- paste("FM", seq(1,5,1), sep = "")
# 

# Scaled value X original sd  + original mean 
# sv * sv.sd + sv.me

### INDIVIDUAL GAMS ####
test <- unique(expand.grid(indicator$INDICATOR, driver$DRIVER[driver$DRIVER != "LANDINGS"],
                           "LANDINGS"))
gams.list <- list()
best.df <- data.frame()
mod.select <- data.frame()
for(ti in 1:nrow(test)) {
  bg.ti <-  test[ti,]
  bg.x <- subset(driver, DRIVER == "LANDINGS")
  bg.y<-  subset(driver, DRIVER == bg.ti[,2])
  bg.z <- subset(indicator, INDICATOR == bg.ti[,1])
  bg.name <- paste(bg.ti[,1],bg.ti[,2], "L", sep = ".")
  bg.year <- intersect(intersect(bg.x$YEAR, bg.y$YEAR), intersect(bg.x$YEAR, bg.z$YEAR))
  bg <- data.frame(x = scale(bg.x$VALUE[bg.x$YEAR %in% bg.year], scale = T, center = T),
                   y = scale(bg.y$VALUE[bg.y$YEAR %in% bg.year], scale = T, center = T),
                   z = scale(bg.z$VALUE[bg.z$YEAR %in% bg.year], scale = T, center = T))
  if(length(bg$x) == length(bg$y)) {
    mat.len <- length(bg$x)
  } else stop("X and Y must be the same length")
  
  # Use GAM to model the relationship between these variables
  k.2 <- ceiling(length(bg.year)/3 - 2)
  k.3 <- ceiling(length(bg.year)/4 - 2)
  
  F1 <- formula(z ~ s(x,  k = k.2))
  F2 <- formula(z ~ s(y,  k = k.2))
  F3 <- formula(z ~ s(x,  k = ceiling(k.2/2)) + s(y,  k = ceiling(k.2/2)))
  F4 <- formula(z ~ s(x,  k = k.2) + s(y,  k = k.2))
  F5 <- formula(z ~ x + s(y,  k = k.2))
  F6 <- formula(z ~ s(x,  k = k.2) + y)
  F7 <- formula(z ~ x + y)
  F8 <- formula(z ~ x)
  F9 <- formula(z ~ y)
  F.list <- paste("F", seq(1,9,1), sep = "")
  
  # Prints a list of all the model output
  gam.ti <- best.gam(formulaList = F.list, dataFrame = bg, gamMethod= "GCV.Cp")
  mod.select <- rbind(mod.select, cbind("ECOSYSTEM" = bg.name, gam.ti))
  
  # Organize models based upon AICc
  best.bic <- gam.ti[order(gam.ti$AICc, decreasing = F),]
 
  
  # Make sure the gam() below matches best.gam(dataFrame, gamMethod)
  best.mod <- gam(get(best.bic[1,1]), method = "GCV.Cp", se = T, data = bg)
  best.df <- rbind(best.df, cbind("INDICATOR_DRIVER" = bg.name, best.bic[1,]))
  gams.list[[bg.name]] <- best.mod
  cat(paste("...", round(ti/nrow(test) * 100), "% ", sep = ""))
}
save(best.df, file = "gamDFATable_v0008.RDATA")
write.csv(best.df, paste0(figure.dir, "dfaTable_v0008.csv"), row.names = F)
# 


# Bivariate GAM smoothers
fin.list <- best.df[best.df$X_EDF > 1 &
                      best.df$Y_EDF > 1 &
                      complete.cases(best.df),]
# fl <- 1
for(fl in 1:nrow(fin.list)) {
  df.name <- as.character(fin.list[fl, 1])
  gam1 <- gams.list[[df.name]]
  
  df.parts <- unlist(strsplit(df.name, "[.]"))
  df.x <- subset(driver, DRIVER == "LANDINGS")
  df.y <- subset(driver, DRIVER == df.parts[2])
  df.z <- subset(indicator, INDICATOR == df.parts[1])
  
  df.year <- intersect(intersect(df.x$YEAR, df.y$YEAR), intersect(df.x$YEAR, df.z$YEAR))
  df.ns <- data.frame(x =  df.x$VALUE[df.x$YEAR %in% df.year],
                      y =  df.y$VALUE[df.y$YEAR %in% df.year],
                      z =  df.z$VALUE[df.z$YEAR %in% df.year])
  
  df <- data.frame(x = scale(df.x$VALUE[df.x$YEAR %in% df.year], scale = T, center = T),
                   y =  scale(df.y$VALUE[df.y$YEAR %in% df.year], scale = T, center = T),
                   z =  scale(df.z$VALUE[df.z$YEAR %in% df.year], scale = T, center = T))
  
  # Create a convex hull surrounding data to aide in plotting
  dia <- max(max(df$x) * .1, max(df$y) * .1)
  hull <- hullFun(df.ns[,c(1:2)], diameter = dia)
  
  if(length(df$x) == length(df$y)) {
    mat.len <- length(df$x)
  } else stop("X and Y must be the same length")
  
  ######################
  # Partial Derivative #
  ######################
  # CONTROL PARAMETERS #
  sp.len <- mat.len * mat.len # Spline length
  # x.ran <- range(df$x)
  # y.ran <- range(df$y)
  x.ran <- range(hull$x)
  y.ran <- range(hull$y)
  eps <- 1e-7 # epsilon for finite difference
  CV <- .95   # critical value for bootstrapped CI
  set.seed(123)
  nb <- 1000  # Bootstrap length
  xseq <- seq(x.ran[1], to = x.ran[2], length.out = mat.len)
  yseq <- seq(y.ran[1], to = y.ran[2], length.out = mat.len)
  eps.list <- matrix(nrow = 9, ncol = 2, 
                     c(+eps, +eps, 0, 0, -eps, 0, -eps, +eps, -eps,
                       +eps, 0, +eps, 0, 0, -eps, -eps, -eps, +eps), 
                     byrow = F)
  k.2 <- ceiling(length(df.year)/3 - 2)
  k.3 <- ceiling(length(df.year)/4 - 2)
  F.i <- summary(gam1)$formula # best.gam() formula
  
  ################### 
  # Bootstrap setup #
  ###################
  #######################
  ### Naive Bootstrap ###
  #######################
  # respboot <- matrix(replicate(nb, df$z[sample(mat.len, rep = TRUE)]), mat.len, nb)
  # drivboot <- matrix(replicate(nb, df$x[sample(mat.len, rep = TRUE)]), mat.len, nb)
  # enviboot <- matrix(replicate(nb, df$y[sample(mat.len, rep = TRUE)]), mat.len, nb)
  
  ################################################
  ### Maximum Entropy Bootstrapped Time Series ###
  ################################################
  respboot <- meboot(df$z, reps = nb, trim = 0.1)$ens 
  drivboot <- meboot(df$x, reps = nb, trim = 0.1)$ens
  enviboot <- meboot(df$y, reps = nb, trim = 0.1)$ens
  
  #############################
  # Empty arrays for new data #
  #############################
  ## Gradient
  grad.Fx <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fy <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fxx <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fyy <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fxy <- array(NA, dim = c(mat.len, mat.len, nb))
  
  grad.Fx.se <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fy.se <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fxx.se <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fyy.se <- array(NA, dim = c(mat.len, mat.len, nb))
  grad.Fxy.se <- array(NA, dim = c(mat.len, mat.len, nb))
  
  ## X and Y range
  dri.ran <- range(drivboot)
  env.ran <- range(enviboot)
  
  ## X and Y sequence
  dri.seq <- seq(dri.ran[1], to = dri.ran[2], length.out = mat.len)
  env.seq <- seq(env.ran[1], to = env.ran[2], length.out = mat.len)
  
  #####################
  # BOOTSTRAP ROUTINE #
  #####################
  for(i in 1:nb) {
    # Isolate the data
    z <- respboot[,i]
    x <- drivboot[,i]
    y <- enviboot[,i]
    # Fit the GAM model
    gam.i <- gam(F.i,  method = "GCV.Cp", se = T)
    # calculate the finite differences for each component of the partial derivative
    tt <- my.grad(the.gam = gam.i, xs = dri.seq, ys = env.seq, eps.mat = eps.list)
    grad.Fx[,,i] <- matrix(tt$x$grad, mat.len, mat.len)  
    grad.Fy[,,i] <- matrix(tt$y$grad, mat.len, mat.len)
    grad.Fxx[,,i] <- matrix(tt$xx$grad, mat.len, mat.len)
    grad.Fyy[,,i] <- matrix(tt$yy$grad, mat.len, mat.len)
    grad.Fxy[,,i] <- matrix(tt$xy$grad, mat.len, mat.len)
    # calculate SE 
    grad.Fx.se[,,i] <- matrix(tt$x$se.grad, mat.len, mat.len)  
    grad.Fy.se[,,i] <- matrix(tt$y$se.grad, mat.len, mat.len)
    grad.Fxx.se[,,i] <- matrix(tt$xx$se.grad, mat.len, mat.len)  
    grad.Fyy.se[,,i] <- matrix(tt$yy$se.grad, mat.len, mat.len)
    grad.Fxy.se[,,i] <- matrix(tt$xy$se.grad, mat.len, mat.len)
    if(i %% 100 == 0) cat(paste("...", i/nb * 100, "% ", sep = ""))
  }
  # 
  
  #######################
  # CALCULATE QUANTILES #
  #######################
  # for each row and column of grad.Fx and grad.Fy take the quantile of nb
  # Set up a matrix for upper
  Fx.CU <- matrix(NA, mat.len, mat.len)
  Fy.CU <- matrix(NA, mat.len, mat.len)
  Fxy.CU <- matrix(NA, mat.len, mat.len)
  Fxx.CU <- matrix(NA, mat.len, mat.len)
  Fyy.CU <- matrix(NA, mat.len, mat.len)
  
  
  # and lower bounds
  Fx.CL <- matrix(NA, mat.len, mat.len)
  Fy.CL <- matrix(NA, mat.len, mat.len)
  Fxy.CL <- matrix(NA, mat.len, mat.len)
  Fxx.CL <- matrix(NA, mat.len, mat.len)
  Fyy.CL <- matrix(NA, mat.len, mat.len)
  
  for(ri in 1:mat.len) {
    for(cj in 1:mat.len) {
      Fx.CU[ri, cj] <- quantile(grad.Fx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
      Fx.CL[ri, cj] <- quantile(grad.Fx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]
      Fy.CU[ri, cj] <- quantile(grad.Fy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
      Fy.CL[ri, cj] <- quantile(grad.Fy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]
      Fxy.CU[ri, cj] <- quantile(grad.Fxy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
      Fxy.CL[ri, cj] <- quantile(grad.Fxy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]
      Fxx.CU[ri, cj] <- quantile(grad.Fxx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
      Fxx.CL[ri, cj] <- quantile(grad.Fxx[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]  
      Fyy.CU[ri, cj] <- quantile(grad.Fyy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[2]
      Fyy.CL[ri, cj] <- quantile(grad.Fyy[ri,cj,], c((1-CV)/2, (1+CV)/2), na.rm = T)[1]  
    }
    cat(paste("[", ri, ",", cj, "] of ", mat.len, "x", mat.len," \n", sep = ""))
  }
  
  
  ####################
  # Predicted values #
  ####################
  
  ## RESPONSE VALUES ##
  respgam <- summary(gam1)$formula # best.gam() formula
  gam2 <- gam(respgam, data = df.ns,   method = "GCV.Cp", se = T)
  
  # xt.ran <- range(df.ns$x)
  # yt.ran <- range(df.ns$y)
  xt.ran <- range(hull$x)
  yt.ran <- range(hull$y)
  
  xtseq <- seq(xt.ran[1], to = xt.ran[2], length.out = mat.len)
  ytseq <- seq(yt.ran[1], to = yt.ran[2], length.out = mat.len)
  
  newDt <- expand.grid(list(x = xtseq, y = ytseq))
  predt <- predict.gam(gam2, newDt, type = "response")
  resp.mat <- matrix(predt, mat.len, mat.len)
  
  ## SCALED AND CENTERED VALUES ##
  pred.grad <- my.grad(the.gam = gam1, xs = xseq, ys = yseq, eps.mat = eps.list)
  newD <- expand.grid(list(x = xseq, y = yseq))
  pred <- predict.gam(gam1, newD, type = "response")
  pred.mat <- matrix(pred, mat.len, mat.len)
  
  Fx.mat <- matrix(pred.grad$x$grad, mat.len, mat.len)
  Fy.mat <- matrix(pred.grad$y$grad, mat.len, mat.len)
  Fxy.mat <- matrix(pred.grad$xy$grad, mat.len, mat.len)
  Fxx.mat <- matrix(pred.grad$xx$grad, mat.len, mat.len)
  Fyy.mat <- matrix(pred.grad$yy$grad, mat.len, mat.len)
  
  est.det <- Fxx.mat * Fyy.mat - (Fxy.mat)^2
  det.CU <- Fxx.CU * Fyy.CU - (Fxy.CU)^2
  det.CL <- Fxx.CL * Fyy.CL - (Fxy.CL)^2
  
  # 
  # hes.mat <- array(NA, c(mat.len, mat.len, 2,2))
  # hes.mat[,,1,1] <- Fxx.mat
  # hes.mat[,,1,2] <- Fxy.mat
  # hes.mat[,,2,1] <- Fxy.mat
  # hes.mat[,,2,2] <- Fyy.mat
  # est.det <- matrix(apply(hes.mat, c(1,2), det), mat.len, mat.len)
  # 
  # hes.CU <- array(NA, c(mat.len, mat.len, 2, 2))
  # hes.CU[,,1,1] <- Fxx.CU
  # hes.CU[,,1,2] <- Fxy.CU
  # hes.CU[,,2,1] <- Fxy.CU
  # hes.CU[,,2,2] <- Fyy.CU
  # 
  # hes.CL <- array(NA, c(mat.len, mat.len, 2, 2))
  # hes.CL[,,1,1] <- Fxx.CL
  # hes.CL[,,1,2] <- Fxy.CL
  # hes.CL[,,2,1] <- Fxy.CL
  # hes.CL[,,2,2] <- Fyy.CL
  # 
  # det.CU <- matrix(apply(hes.CU, c(1,2), det), mat.len, mat.len)
  # det.CL <- matrix(apply(hes.CL, c(1,2), det), mat.len, mat.len)
  
  save(gam1,                          # Original model
       df, df.ns,                     # S&C and Original data      
       xtseq, ytseq,                  # Original seq
       xseq, yseq,                    # S&C seq   
       grad.Fx, grad.Fy, grad.Fxx, grad.Fyy, grad.Fxy, #
       Fx.CU, Fx.CL,                  #  
       Fy.CU, Fy.CL,                  #
       Fxx.CU,Fxx.CL,                 #    
       Fyy.CU,Fyy.CL,                 # 
       Fxy.CU,Fxy.CL,                 #
       df.parts, df.name,             # Name controls
       pred.mat,                      # S&C Predicted matrix 
       resp.mat,                      # Response predicted matrix 
       det.CU, det.CL,                # Determinant CI
       hull,                          # Convex hull
       file = paste0(df.name, "_DATA_v0008.RDATA"))
}

###########
## PLOTS ##
###########
col.mat <- matrix(NA, mat.len, mat.len)
col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CL > 0] <- 1 # Local minimum
col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CU < 0] <- 2 # Local maximum
col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CU < 0] <- 3 # Saddle point

grad.plot <- paste(figure.dir, df.name, "v_002.png", sep="")
png(file = grad.plot, width = 76, height = 80, units = "mm", res = 600)
# layout(matrix(c(1,2), 1,2,  byrow = T),
#        widths = lcm(c(8.45, 8.45)), heights = lcm(c(4.75, 4.75)),
#        respect = F)
nlev <-  10
levs <- pretty(range(pred.mat), nlev)

par(mar=c(2.15,2.15,1.5,0.25), 
    #     oma = c(0,0,0,0),
    mgp=c(1,0.5,0), 
    cex= .75, 
    tck= -0.015, 
    family = "sans")
plot(df$x, df$y, type = "n",
     ylab = "",
     xlab = "")
image(xseq,
      yseq,
      pred.mat,
      col = gray.colors(50, start = 0.5, end = 1, gamma = 1, alpha = NULL),
      add = T)
image(xseq,
      yseq,
      col.mat,
      #       col = "white",
      col= c("gray20", "gray30", "gray40"), 
      add = T)
# contour(x = xseq,
#         y = yseq,
#         z = pred.mat,
#         add = T, drawlabels = F)

contour(x = xseq,
        y = yseq,
        z = pred.mat,
        zlim = range(pred.mat, finite = T),
        lwd = seq(0, 2, length.out= length(levs)),
        #         col = gray(10:0/ 11),
        col = "black",
        drawlabels = T,
        method = "edge",
        #         labcex = 1,
        #         vfont = c("serif", "bold"),
        vfont = NULL,
        add = T)
legend("topleft", 
       legend = df.parts[1],
       bty = "n")
dev.off()

# 
# col.mat <- matrix(NA, mat.len, mat.len)
# col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CL > 0] <- "red" # Local minimum
# col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CU < 0] <-  "black" # Local maximum
# col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CU < 0] <- "yellow" # Saddle point
# 
# open3d()
# persp3d(xseq, 
#         yseq, 
#         pred.mat, 
#         col = col.mat, 
#         xlab = "LANDINGS", 
#         ylab = df.parts[2],
#         zlab = df.parts[1])
# rgl.spheres(df$x,
#             df$y, 
#             df$z, 
#             color = "black",
#             radius = .05)
# par3d(params = list(
#   windowRect = c(100, 100, 600, 600)))
# movie3d(spin3d(axis = c(0,0,1),
#                rpm=10),
#         duration=6, 
#         movie = paste0("movie", df.name), 
#         dir = figure.dir,
#         clean = TRUE)
# rgl.close()
# }




# rbPal <- colorRampPalette(c('white','black'))
# df$Col <- rbPal(10)[as.numeric(cut(df$z,breaks = 10))]


# # 
# # col.mat <- matrix(NA, mat.len, mat.len)
# # col.mat[Fx.CL > 0] <- "red4" # fx as y is held constant (pos trend)
# # col.mat[Fy.CL > 0] <- "red1" # fy as x is held constant (pos trend)
# # col.mat[Fx.CU < 0] <- "green4" # fx as y is held constant (neg trend)
# # col.mat[Fy.CU < 0] <- "green1" # fy as s is held constant (neg trend)
# # 
# col.mat <- matrix(NA, mat.len, mat.len)
# col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CL > 0] <- "gray80" # Local minimum
# col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CL > 0 & Fxx.CU < 0] <- "gray60" # Local maximum
# col.mat[Fx.CL < 0 & Fy.CL < 0 & Fx.CU > 0 & Fy.CU > 0 & det.CU < 0] <- "gray40" # Saddle point
# 
# # perspPar <- par3d(no.readonly=TRUE)
# # save(perspPar, file = "perspPlotPar.Rdata")
# load("perspPlotPar.Rdata") # sets the rotation parameters
# # perspPar$family <- "serif"
# # perspPar$font <- 1
# # perspPar$useFreeType = F
# par(mfg= c(1,2))
# perspPlot <- persp3d(xseq, 
#                      yseq, 
#                      pred.mat, 
#                      col = col.mat, 
#                      xlab = "Landings",
#                      ylab = df.parts[2],
#                      zlab = df.parts[1],
#                      box = F)
# par3d(perspPar)
# # rgl.postscript(paste(df.name, "pdf", sep = "."),"pdf",drawText=TRUE)
# # rgl.snapshot(paste(df.name, "png", sep = "."), "png")
# rgl.close()
# # par3d()$family


########################################################################################
# If f(x,y)  is a two-dimensional function that has a relative extremum at a point  
# (x0, y0) and has continuous partial derivatives at this point, 
# then fx(x0,y0) = 0 and fy(x0,y0) = 0. The second partial derivatives test classifies 
# the point as a local maximum or relative minimum.

## For the determinant M(x,y) = det(H(x,y)) = 
# if (a, b) is a critical point of f (that is, fx(a, b) = fy(a, b) = 0)
# i.e., fx.CL < 0 & fy.CL < 0 & fx.CU > 0 & fy.CU >0 # We have 95% certainty that we contain zero...

# 1. if M(a,b) > 0 and Fxx(ab) > 0 then (a,b) is a local minimum of f
# 2. if M(a,b) > 0 and Fxx(ab) < 0 then (a,b) is a local maximum of f
# 3. if M(a,b) < 0 and Fxx(ab) then (a,b) is a saddle point of f
# 4. if M(a,b) = 0 then the second derivative test is inconclusive, 
#    and the point (a,b) could be any of a minimum, maximum or saddle point
# Fxx | Fxy
# Fyx | Fyy


persp3d(xseq, yseq, Fx.mat, col = col.mat)
persp3d(xseq, yseq, Fx.CL, col = "RED", add = F)
persp3d(xseq, yseq, Fx.CU, col = "GREEN", add = T)
planes3d(0,0,1)


persp3d(xseq, yseq, Fy.mat, col = col.mat)
persp3d(xseq, yseq, Fy.CL, col = "RED", add = T)
persp3d(xseq, yseq, Fy.CU, col = "GREEN", add = T)
planes3d(0,0,1)

persp3d(xseq, yseq, est.det, col = col.mat)
persp3d(xseq, yseq, det.CL, col = "RED", add = F)
persp3d(xseq, yseq, det.CU, col = "GREEN", add = T)
planes3d(0,0,1)


