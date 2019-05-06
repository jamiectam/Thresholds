################
# FORMULA INFO #
################
# PURPOSE: Function that uses BIC model selection to identify the "best" model for a data series
# and creates a data.frame with with all relevant model outputs.
# AUTHOR: Scott Large 2013
# REVIEWED BY:
# VERSION: 0.1
#

# formulaList is a list(c("F1"...)) of GAM models (where, F1 <- formula(z ~ s(x,  k = k.2))...)
# dataFrame is the data.frame to be analyzed, where df <- data.frame("x" = x, y = y, "z" = z)
# gamMethod (see ?mgcv for description of methods) "GCV.Cp", "GACV.Cp", "REML", "P-REML", "ML", and "P-ML"

best.gam <- function(formulaList, dataFrame, gamMethod) {
  # Make sure "mgcv" is loaded:
  if(require("mgcv", character.only = T)) {
  } else {
    print("trying to install mgcv")
    install.packages("mgcv")
    if(require(mgcv)) {
      print("mgcv installed and loaded")
    } else {
      stop("could not install mgcv")
    }
  }
  gam.fits <- data.frame()
  for(f in formulaList){
    gam.f <- gam(get(f), method = gamMethod, se = T, data = dataFrame)
    ## Collect appropriate P-values and estimated degrees of freedom for each term
    # For models with 2 smoothing terms
    if(length(summary(gam.f)$s.pv) == 2) { 
      X_pvalue <- round(summary(gam.f)$s.pv, 3)[1]
      Y_pvalue <- round(summary(gam.f)$s.pv, 3)[2]
      X_EDF    <- round(summary(gam.f)$edf, 3)[1]
      Y_EDF    <- round(summary(gam.f)$edf, 3)[2]
    } 
    # For models with linear and nonlinear components
    if(length(summary(gam.f)$s.pv) == 1 & length(summary(gam.f)$pTerms.pv) == 1) { 
      # Linear "x" and nonlinear "y"
      if(all(dimnames(summary(gam.f)$pTerms.pv) == "x")) {
        X_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
        Y_pvalue <- round(summary(gam.f)$s.pv, 3)
        X_EDF    <- NA
        Y_EDF    <- round(summary(gam.f)$edf, 3)
      }
      # Nonlinear "x" and linear "y"
      if(all(dimnames(summary(gam.f)$pTerms.pv) == "y")) {
        X_pvalue <- round(summary(gam.f)$s.pv, 3)
        Y_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
        X_EDF    <- round(summary(gam.f)$edf, 3)
        Y_EDF    <- NA
      } 
    }
    # For models with linear "x" and "y"
    if(length(summary(gam.f)$pTerms.pv) == 2) {
      X_pvalue <- round(summary(gam.f)$pTerms.pv, 3)[1]
      Y_pvalue <- round(summary(gam.f)$pTerms.pv, 3)[2]
      X_EDF    <- NA
      Y_EDF    <- NA
    }
    # For models with a single driver: s(x), s(y), x, and y
    if(length(summary(gam.f)$s.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "x")) {
      X_pvalue <- round(summary(gam.f)$s.pv, 3)
      Y_pvalue <- NA
      X_EDF    <- round(summary(gam.f)$edf, 3)
      Y_EDF    <- NA
    }
    if(length(summary(gam.f)$s.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "y")) {
      X_pvalue <- NA
      Y_pvalue <- round(summary(gam.f)$s.pv, 3)
      X_EDF    <- NA
      Y_EDF    <- round(summary(gam.f)$edf, 3)
    }
    if(length(summary(gam.f)$pTerms.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "x")) {
      X_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
      Y_pvalue <- NA
      X_EDF    <- NA
      Y_EDF    <- NA
    }
    if(length(summary(gam.f)$pTerms.pv) == 1 & all(attr(terms(gam.f), "term.labels") == "y")) {
      X_pvalue <- NA
      Y_pvalue <- round(summary(gam.f)$pTerms.pv, 3)
      X_EDF    <- NA
      Y_EDF    <- NA
    }
    ## Create a data frame of all GAM for formulaList
    gam.fits <- rbind(gam.fits,
                      data.frame("MODEL" = f,
                                 "AICc"= AICc(gam.f),
                                 "BIC" = BIC(gam.f),
                                 "GCF" = summary(gam.f)$sp.criterion,
                                 "DEV" = summary(gam.f)$dev.expl,
                                 "X_p-value" = X_pvalue, 
                                 "Y_p-value" = Y_pvalue,
                                 "X_EDF" = X_EDF,
                                 "Y_EDF" = Y_EDF,
                                 stringsAsFactors=FALSE))
    assign(paste(f, "gam", sep = "."), gam.f)
  }
  return(gam.fits)
}

# Calculate AICc, which corrects for finite sample size Burnham & Anderson 2002
AICc <- function(mod) {
  K.c <- mod$rank
  N.c <- length(mod$residuals)
  AIC.c <- round(mod$aic + (2*K.c*(K.c+1)/(N.c-K.c-1)),3)
  return(AIC.c)
}

#################
# FUNCTION INFO #
#################
# PURPOSE: Create a convex hull around raw data with a buffer of 1/diameter 
# AUTHOR: Scott Large 2013
# REVIEWED BY:
# vERSION: 0.1
#

hullFun <- function(hull.df, diameter = 1, npoints = 10) {
  
  circleFun <- function(center = c(0,0), diameter = 1, npoints = 10){
    r = diameter / 2
    tt <- seq(0,2*pi,length.out = npoints)
    xx <- center[1] + r * cos(tt)
    yy <- center[2] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
  }
  
  # Scale and center original values
  scaled.df <- data.frame(x = scale(hull.df[,1], scale = T, center = T),
                          y = scale(hull.df[,2], scale = T, center = T))
  
  # Take convex hull of the scaled and centered values
  hpts <- chull(scaled.df)
  hpts <- c(hpts, hpts[1])
  # Create a circle of points surrounding each convex hull point
  circ <- data.frame()
  for(cf in 1:length(hpts)) {
    cf.df <- c(scaled.df[hpts[cf],1], scaled.df[hpts[cf],2])
    circ <- rbind(circ, circleFun(center = cf.df, diameter = diameter, npoints = npoints))
  }
  
  hpts2 <- chull(circ)
  hpts2 <- c(hpts2,hpts2[1])
  hull <- data.frame(x = circ[hpts2,1] * sd(hull.df[,1]) + mean(hull.df[,1]),
                     y = circ[hpts2,2] * sd(hull.df[,2]) + mean(hull.df[,2]))
  return(hull)
}

#################
# FUNCTION INFO #
#################
# PURPOSE: Uses the predicted values of "the.gam" to calculate partial derivatives of x and y (finite differences method)
# AUTHOR: Scott Large 2013
# REVIEWED BY:
# vERSION: 0.1
#

# the.gam <- gam1
my.grad <- function(the.gam, xs, ys, eps.mat) { 
  # Creates an array of predicted values from "the.gam" along x and y sequences "xs", and "ys", for
  # each row in the esp.mat. For partial derivatives of x and y use eps.mat listed below: # 
  # eps <- 1e-7
  # eps.list <- matrix(nrow = 9, ncol = 2, c(+eps, +eps, 0, 0, -eps, 0, -eps, +eps, -eps,
  #                                          +eps, 0, +eps, 0, 0, -eps, -eps, -eps, +eps), byrow = F)
  ### CREATE EPS TERMS ###  
  m.terms <- attr(terms(the.gam), "term.labels") # Select model terms i.e. "X" and "Y" or "Env" and "Dri"
  lp.terms <- the.gam$coefficients # coefficients in the model, to create an array
  p.array <- array(NA, dim = c(sp.len, length(lp.terms), nrow(eps.list))) # sp.len X #coeffecients X eps terms array
  for(i in 1:nrow(eps.mat)) {
    new.dd <- expand.grid(list(x = xs + eps.mat[i,1], # Creates an expanded matrix based on the x and y sequences
                               y = ys + eps.mat[i,2]))
    colnames(new.dd) <- m.terms # Adds the appropriate column names for the gam model
    p.i <- predict.gam(the.gam, newdata = new.dd, type = "lpmatrix") # predicts using the gam for the provided data
    p.array[,,i] <- p.i  # Puts each eps term in the array
  }
  # Create list to hold the data
  nt <- length(m.terms)
  lD <- vector(mode = "list", length = nt + 3)
  names(lD) <- c(m.terms, "xx", "yy", "xy")
  ### CALCULATE Fx, Fx, Fxx, Fyy, and Fxy ###
  Xp <- (p.array[,,2] - p.array[,,5]) / (2*eps) # f(x+h, y) - f(x-h, y)/2h
  Yp <- (p.array[,,3] - p.array[,,6]) / (2*eps) # f(x, y+k) - f(x, y-k)/2k
  XXp <-(p.array[,,2] - 2*(p.array[,,4]) +  p.array[,,5]) / (eps^2) # f(x+h, y) - 2f(x,y) + f(x-h, y) / h^2
  YYp <-(p.array[,,3] - 2*(p.array[,,4]) +  p.array[,,6]) / (eps^2) # f(x, y+h) - 2f(x,y) + f(x, y-h) / k^2
  # f(x+h, y+k) - f(x+h, y) - f(x, y+k) + 2f(x,y) - f(x-h, y) - f(x, y-k) + f(x-h, y-k))/(2hk)
  XYp <- (p.array[,,1] - p.array[,,2] -p.array[,,3] + (2*p.array[,,4]) - p.array[,,5] - p.array[,,6] + p.array[,,7])/ (2*eps*eps) 
  
  ### MATRIX MULTIPLY THE COEF ###
  Xi <- Xp * 0
  x.want <- grep(m.terms[1], names(lp.terms))
  Xi[, x.want] <- Xp[, x.want]
  x.df <- Xi %*% coef(the.gam)
  x.df.sd <- rowSums(Xi %*% the.gam$Vp * Xi)^.5
  lD[[1]] <- list(grad = x.df, se.grad = x.df.sd)
  # and for Y's
  Yi <- Yp * 0
  y.want <- grep(m.terms[2], names(lp.terms))
  Yi[, y.want] <- Yp[, y.want]
  y.df <- Yi %*% coef(the.gam)
  y.df.sd <- rowSums(Yi %*% the.gam$Vp * Yi)^.5
  lD[[2]] <- list(grad = y.df, se.grad = y.df.sd)
  # # and for Fxx
  XXi <- XXp * 0
  xx.want <- grep(m.terms[1], names(lp.terms))
  XXi[, xx.want] <- XXp[, xx.want]
  xx.df <- XXi %*% coef(the.gam)
  xx.df.sd <- rowSums(XXi %*% the.gam$Vp * XXi)^.5
  lD[[3]] <- list(grad = xx.df, se.grad = xx.df.sd)
  # # and for Fyy
  YYi <- YYp * 0
  yy.want <- grep(m.terms[2], names(lp.terms))
  YYi[, yy.want] <- YYp[, yy.want]
  yy.df <- YYi %*% coef(the.gam)
  yy.df.sd <- rowSums(YYi %*% the.gam$Vp * YYi)^.5
  lD[[4]] <- list(grad = yy.df, se.grad = yy.df.sd)
  # # and for Fxy
  XYi <- XYp * 0
  XYi[, -1] <- XYp[,-1]
  xy.df <- XYi %*% coef(the.gam)
  xy.df.sd <- rowSums(XYi %*% the.gam$Vp * XYi)^.5
  lD[[5]] <- list(grad = xy.df, se.grad = xy.df.sd)
  # Package it all up
  class(lD) <- "my.grad"
  #   lD$gamModel <- the.gam
  #   lD$eps <- eps
  #   lD$eval <- expand.grid(list(x = xs, y = ys))
  return(lD)
}
