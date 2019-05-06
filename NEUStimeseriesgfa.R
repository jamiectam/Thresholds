###Time series experiment### Is 10 years enough to capture tresholds?
##Is the longterm data caputuring the same threshold points as short term?
##threshold changes in ecosystem status are more useful than baselines. 
##Using shifts in pressure response relationships are more useful than historical baselines


set.seed(627)
#########
require(reshape2)
require(gradientForest)
require(ggplot2)
require(extendedForest)
#################
#create dataframes for NEUSenv and NEUSbio
NEUS10bio <- NEUSbio[40:50,]
NEUS20bio <- NEUSbio[30:50,]
NEUS30bio <- NEUSbio[20:50,]
NEUS40bio <- NEUSbio[10:50,]
NEUS50bio <- NEUSbio

NEUS10env <- NEUSenv[40:50,]
NEUS20env <- NEUSenv[30:50,]
NEUS30env <- NEUSenv[20:50,]
NEUS40env <- NEUSenv[10:50,]
NEUS50env <- NEUSenv



##########Gradient Forest
#Set up data (use region specific env and bio dataframes)
dat10<-cbind(NEUS10env, NEUS10bio)
dat20<-cbind(NEUS20env, NEUS20bio)
dat30<-cbind(NEUS30env, NEUS30bio)
dat40<-cbind(NEUS40env, NEUS40bio)
dat50<-cbind(NEUS50env, NEUS50bio)

ind.name10<-colnames(NEUS10bio)
ind.name20<-colnames(NEUS20bio)
ind.name30<-colnames(NEUS30bio)
ind.name40<-colnames(NEUS40bio)
ind.name50<-colnames(NEUS50bio)

dri.name10<-colnames(NEUS10env)
dri.name20<-colnames(NEUS20env)
dri.name30<-colnames(NEUS30env)
dri.name40<-colnames(NEUS40env)
dri.name50<-colnames(NEUS50env)


# Many of the pressure variables do not have full data (some NA):
## Option 1- Impute:
dat10 <- na.roughfix(dat10)
dat20 <- na.roughfix(dat20)
dat30 <- na.roughfix(dat30)
dat40 <- na.roughfix(dat40)
dat50 <- na.roughfix(dat50)

#
## Option 2- Get rid of the columns without full time series:
#
# dat.sc <- dat.full[,apply(dat.full, 2, function(x)!all(is.na(x)))]
#
# Somewhat roundabout way of doing this, but can also be used to find the best
# contiguous set of data.
# if any column has less than XX years, omit.
#
# cuts <- 10
# len.list <- sapply(dat.sc, function(x) length(na.contiguous(x)))
# keep.list <- names(len.list[len.list >= cuts])
# dat.kl <- dat.sc[, keep.list] 
# dat <- dat.kl[apply(dat.kl, 1, function(x)!any(is.na(x))),]
#
## Option 3- Simply use as is (what I initially tried, but it might not run...)
# dat <- dat.full
#dat1<-dat.full



# Maximum level of splits
lev10 <- floor(log2(nrow(dat10) * 0.368/2))
#colnames(dat)[!colnames(dat) %in% c(dri.name, ind.name)]\
#dat<-data.frame(dat)
#
## GF analysis ##
NEUS10gfa <-   gradientForest(data = dat10, 
                            predictor.vars = dri.name10, 
                            response.vars = ind.name10,
                            ntree = 1000, 
                            transform = NULL,
                            maxLevel = 1,
                            corr.threshold = 0.5, 
                            compact = F,
                            trace = T)

#save Rdata

save(NEUS10gfa, file = "NEUS10gfa_v01.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(NEUS10gfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- NEUS10gfa$imp.rsq
indicatorNA <- colnames(indicatorIMP)[apply(indicatorIMP, 2, function(x)all(is.na(x)))]
indicatorIMP <- indicatorIMP[apply(indicatorIMP, 2, function(x)!all(is.na(x))),]
#
indicatorDF <- data.frame("VARIABLE" = colnames(indicatorIMP),
                          "TYPE" = "INDICATOR", melt(apply(indicatorIMP, 2, mean, na.rm = T),
                                                     variable.name = "VARIABLE",
                                                     value.name = "MEAN"),
                          melt(apply(indicatorIMP, 2, min, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MIN"),
                          melt(apply(indicatorIMP, 2, max, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MAX"))
#
pressureIMP <- NEUS10gfa$imp.rsq
pressureNA <- row.names(pressureIMP)[apply(pressureIMP, 1, function(x)all(is.na(x)))]
pressureIMP <- pressureIMP[apply(pressureIMP, 1, function(x)!all(is.na(x))),]
pressureDF <- data.frame("VARIABLE" = row.names(pressureIMP),
                         "TYPE" = "PRESSURE", melt(apply(pressureIMP, 1, mean, na.rm = T),
                                                   variable.name = "VARIABLE",
                                                   value.name = "MEAN"),
                         melt(apply(pressureIMP, 1, min, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MIN"),
                         melt(apply(pressureIMP, 1, max, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MAX"))
#
indicatorDF[,c(3:5)] <- apply(indicatorDF[,c(3:5)], 2, signif, 2)
pressureDF[,c(3:5)] <- apply(pressureDF[,c(3:5)], 2, signif, 2)
# indicatorDF$VARIABLE <- factor(indicatorDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))
# pressureDF$VARIABLE <- factor(pressureDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))

#pressureDF <- rbind(pressureDF, cbind(VARIABLE = pressureNA,TYPE = "PRESSURE", MEAN = "NA", MIN = "NA", MAX = "NA"))
#
modPerformance <- rbind(pressureDF, indicatorDF)
##
write.csv(modPerformance, file = "NEUS10gfModelPerformance_v01.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(NEUS10gfa)[importance(NEUS10gfa) > 0])

#
# Overall Importance
png(file = "NEUS10Importance_v01.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(NEUS10gfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "NEUS10splitRatio_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(NEUS10gfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "NEUS10splitImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(NEUS10gfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "NEUS10splitData_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(NEUS10gfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="NEUS10cumImportanceSplit_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS10gfa, 
          imp.vars= imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 6,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = "NEUS10cumImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS10gfa, 
          imp.vars = imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 6,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(NEUS10gfa, dat[, imp.vars])
row.names(Trns_grid) <- c(2004:2013)
NEUS10PCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="NEUS10PCA_v01.png", PCbiplot(NEUS10PCs),
       width = 83, height = 83, units = "mm", scale = 2)
###########

#find names of important variables
imp.vars
#pull data to get threshold quantiles from split density
a<-NEUS10gfa$d$Exploitation$x
b<-NEUS10gfa$d$Population$x
c<-NEUS10gfa$d$Landings$x
d<-NEUS10gfa$d$AMO$x
e<-NEUS10gfa$d$SST$x
f<-NEUS10gfa$d$Seafood$x
g<-NEUS10gfa$d$NAO_w$x
h<-NEUS10gfa$d$GS$x
i<-NEUS10gfa$d$MEI$x
NEUS.dens10<-cbind(a,b,c,d,e,f,g,h,i)

NEUS.dens10<-as.data.frame(NEUS.dens10)
colnames(NEUS.dens10)<-c("Exploitation", "Population", "Landings", "AMO", "SST", "Seafood", "NAO_w", "GS", "MEI")
neus.d10<-sapply(NEUS.dens10, summary)
write.csv(neus.d10, file="NEUS.dens10.split.csv")


##############
#20 years ####
##############

lev20 <- floor(log2(nrow(dat20) * 0.368/2))
#colnames(dat)[!colnames(dat) %in% c(dri.name, ind.name)]\
#dat<-data.frame(dat)
#
## GF analysis ##
NEUS20gfa <-   gradientForest(data = dat20, 
                              predictor.vars = dri.name20, 
                              response.vars = ind.name20,
                              ntree = 1000, 
                              transform = NULL,
                              maxLevel = lev20,
                              corr.threshold = 0.5, 
                              compact = F,
                              trace = T)

#save Rdata

save(NEUS20gfa, file = "NEUS20gfa_v01.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(NEUS20gfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- NEUS20gfa$imp.rsq
indicatorNA <- colnames(indicatorIMP)[apply(indicatorIMP, 2, function(x)all(is.na(x)))]
indicatorIMP <- indicatorIMP[apply(indicatorIMP, 2, function(x)!all(is.na(x))),]
#
indicatorDF <- data.frame("VARIABLE" = colnames(indicatorIMP),
                          "TYPE" = "INDICATOR", melt(apply(indicatorIMP, 2, mean, na.rm = T),
                                                     variable.name = "VARIABLE",
                                                     value.name = "MEAN"),
                          melt(apply(indicatorIMP, 2, min, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MIN"),
                          melt(apply(indicatorIMP, 2, max, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MAX"))
#
pressureIMP <- NEUS20gfa$imp.rsq
pressureNA <- row.names(pressureIMP)[apply(pressureIMP, 1, function(x)all(is.na(x)))]
pressureIMP <- pressureIMP[apply(pressureIMP, 1, function(x)!all(is.na(x))),]
pressureDF <- data.frame("VARIABLE" = row.names(pressureIMP),
                         "TYPE" = "PRESSURE", melt(apply(pressureIMP, 1, mean, na.rm = T),
                                                   variable.name = "VARIABLE",
                                                   value.name = "MEAN"),
                         melt(apply(pressureIMP, 1, min, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MIN"),
                         melt(apply(pressureIMP, 1, max, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MAX"))
#
indicatorDF[,c(3:5)] <- apply(indicatorDF[,c(3:5)], 2, signif, 2)
pressureDF[,c(3:5)] <- apply(pressureDF[,c(3:5)], 2, signif, 2)
# indicatorDF$VARIABLE <- factor(indicatorDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))
# pressureDF$VARIABLE <- factor(pressureDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))

#pressureDF <- rbind(pressureDF, cbind(VARIABLE = pressureNA,TYPE = "PRESSURE", MEAN = "NA", MIN = "NA", MAX = "NA"))
#
modPerformance <- rbind(pressureDF, indicatorDF)
##
write.csv(modPerformance, file = "NEUS20gfModelPerformance_v01.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(NEUS20gfa)[importance(NEUS20gfa) > 0])

#
# Overall Importance
png(file = "NEUS20Importance_v01.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(NEUS20gfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "NEUS20splitRatio_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(NEUS20gfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "NEUS20splitImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(NEUS20gfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "NEUS20splitData_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(NEUS20gfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="NEUS20cumImportanceSplit_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS20gfa, 
          imp.vars= imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 6,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = "NEUS20cumImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS20gfa, 
          imp.vars = imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 6,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(NEUS20gfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1994:2013)
NEUS20PCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="NEUS20PCA_v01.png", PCbiplot(NEUS20PCs),
       width = 83, height = 83, units = "mm", scale = 2)
###########

#find names of important variables
imp.vars
#pull data to get threshold quantiles from split density
a<-NEUS20gfa$d$Exploitation$x
b<-NEUS20gfa$d$Population$x
c<-NEUS20gfa$d$Landings$x
d<-NEUS20gfa$d$AMO$x
e<-NEUS20gfa$d$SST$x
f<-NEUS20gfa$d$Seafood$x
g<-NEUS20gfa$d$NAO_w$x
h<-NEUS20gfa$d$GS$x
i<-NEUS20gfa$d$MEI$x
NEUS20.dens<-cbind(a,b,c,d,e,f,g,h,i)

NEUS20.dens<-as.data.frame(NEUS20.dens)
colnames(NEUS20.dens)<-c("Exploitation", "Population", "Landings", "AMO", "SST", "Seafood", "NAO_w", "GS", "MEI")
neus20.d<-sapply(NEUS20.dens, summary)
write.csv(neus20.d, file="NEUS20.dens.split.csv")

######################
####30years#########
###################

lev30 <- floor(log2(nrow(dat30) * 0.368/2))
#colnames(dat)[!colnames(dat) %in% c(dri.name, ind.name)]\
#dat<-data.frame(dat)
#
## GF analysis ##
NEUS30gfa <-   gradientForest(data = dat30, 
                              predictor.vars = dri.name30, 
                              response.vars = ind.name30,
                              ntree = 1000, 
                              transform = NULL,
                              maxLevel = lev30,
                              corr.threshold = 0.5, 
                              compact = F,
                              trace = T)

#save Rdata

save(NEUS30gfa, file = "NEUS30gfa_v01.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(NEUS30gfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- NEUS30gfa$imp.rsq
indicatorNA <- colnames(indicatorIMP)[apply(indicatorIMP, 2, function(x)all(is.na(x)))]
indicatorIMP <- indicatorIMP[apply(indicatorIMP, 2, function(x)!all(is.na(x))),]
#
indicatorDF <- data.frame("VARIABLE" = colnames(indicatorIMP),
                          "TYPE" = "INDICATOR", melt(apply(indicatorIMP, 2, mean, na.rm = T),
                                                     variable.name = "VARIABLE",
                                                     value.name = "MEAN"),
                          melt(apply(indicatorIMP, 2, min, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MIN"),
                          melt(apply(indicatorIMP, 2, max, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MAX"))
#
pressureIMP <- NEUS30gfa$imp.rsq
pressureNA <- row.names(pressureIMP)[apply(pressureIMP, 1, function(x)all(is.na(x)))]
pressureIMP <- pressureIMP[apply(pressureIMP, 1, function(x)!all(is.na(x))),]
pressureDF <- data.frame("VARIABLE" = row.names(pressureIMP),
                         "TYPE" = "PRESSURE", melt(apply(pressureIMP, 1, mean, na.rm = T),
                                                   variable.name = "VARIABLE",
                                                   value.name = "MEAN"),
                         melt(apply(pressureIMP, 1, min, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MIN"),
                         melt(apply(pressureIMP, 1, max, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MAX"))
#
indicatorDF[,c(3:5)] <- apply(indicatorDF[,c(3:5)], 2, signif, 2)
pressureDF[,c(3:5)] <- apply(pressureDF[,c(3:5)], 2, signif, 2)
# indicatorDF$VARIABLE <- factor(indicatorDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))
# pressureDF$VARIABLE <- factor(pressureDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))

#pressureDF <- rbind(pressureDF, cbind(VARIABLE = pressureNA,TYPE = "PRESSURE", MEAN = "NA", MIN = "NA", MAX = "NA"))
#
modPerformance <- rbind(pressureDF, indicatorDF)
##
write.csv(modPerformance, file = "NEUS30gfModelPerformance_v01.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(NEUS30gfa)[importance(NEUS30gfa) > 0])

#
# Overall Importance
png(file = "NEUS30Importance_v01.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(NEUS30gfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "NEUS30splitRatio_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(NEUS30gfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "NEUS30splitImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(NEUS30gfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "NEUS30splitData_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(NEUS30gfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="NEUS30cumImportanceSplit_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS30gfa, 
          imp.vars= imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 6,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = "NEUS30cumImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS30gfa, 
          imp.vars = imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 6,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(NEUS30gfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1984:2013)
NEUS30PCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="NEUS30PCA_v01.png", PCbiplot(NEUS30PCs),
       width = 83, height = 83, units = "mm", scale = 2)
###########

#find names of important variables
imp.vars
#pull data to get threshold quantiles from split density
a<-NEUS30gfa$d$Exploitation$x
b<-NEUS30gfa$d$Population$x
c<-NEUS30gfa$d$Landings$x
d<-NEUS30gfa$d$AMO$x
e<-NEUS30gfa$d$SST$x
f<-NEUS30gfa$d$Seafood$x
g<-NEUS30gfa$d$NAO_w$x
h<-NEUS30gfa$d$GS$x
i<-NEUS30gfa$d$MEI$x
NEUS30.dens<-cbind(a,b,c,d,e,f,g,h,i)

NEUS30.dens<-as.data.frame(NEUS30.dens)
colnames(NEUS30.dens)<-c("Exploitation", "Population", "Landings", "AMO", "SST", "Seafood", "NAO_w", "GS", "MEI")
NEUS30.d<-sapply(NEUS30.dens, summary)
write.csv(NEUS30.d, file="NEUS30.dens.split.csv")

###################
###40 years######
####################


lev40 <- floor(log2(nrow(dat40) * 0.368/2))
#colnames(dat)[!colnames(dat) %in% c(dri.name, ind.name)]\
#dat<-data.frame(dat)
#
## GF analysis ##
NEUS40gfa <-   gradientForest(data = dat40, 
                              predictor.vars = dri.name40, 
                              response.vars = ind.name40,
                              ntree = 1000, 
                              transform = NULL,
                              maxLevel = lev40,
                              corr.threshold = 0.5, 
                              compact = F,
                              trace = T)

#save Rdata

save(NEUS40gfa, file = "NEUS40gfa_v01.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(NEUS40gfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- NEUS40gfa$imp.rsq
indicatorNA <- colnames(indicatorIMP)[apply(indicatorIMP, 2, function(x)all(is.na(x)))]
indicatorIMP <- indicatorIMP[apply(indicatorIMP, 2, function(x)!all(is.na(x))),]
#
indicatorDF <- data.frame("VARIABLE" = colnames(indicatorIMP),
                          "TYPE" = "INDICATOR", melt(apply(indicatorIMP, 2, mean, na.rm = T),
                                                     variable.name = "VARIABLE",
                                                     value.name = "MEAN"),
                          melt(apply(indicatorIMP, 2, min, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MIN"),
                          melt(apply(indicatorIMP, 2, max, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MAX"))
#
pressureIMP <- NEUS40gfa$imp.rsq
pressureNA <- row.names(pressureIMP)[apply(pressureIMP, 1, function(x)all(is.na(x)))]
pressureIMP <- pressureIMP[apply(pressureIMP, 1, function(x)!all(is.na(x))),]
pressureDF <- data.frame("VARIABLE" = row.names(pressureIMP),
                         "TYPE" = "PRESSURE", melt(apply(pressureIMP, 1, mean, na.rm = T),
                                                   variable.name = "VARIABLE",
                                                   value.name = "MEAN"),
                         melt(apply(pressureIMP, 1, min, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MIN"),
                         melt(apply(pressureIMP, 1, max, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MAX"))
#
indicatorDF[,c(3:5)] <- apply(indicatorDF[,c(3:5)], 2, signif, 2)
pressureDF[,c(3:5)] <- apply(pressureDF[,c(3:5)], 2, signif, 2)
# indicatorDF$VARIABLE <- factor(indicatorDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))
# pressureDF$VARIABLE <- factor(pressureDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))

#pressureDF <- rbind(pressureDF, cbind(VARIABLE = pressureNA,TYPE = "PRESSURE", MEAN = "NA", MIN = "NA", MAX = "NA"))
#
modPerformance <- rbind(pressureDF, indicatorDF)
##
write.csv(modPerformance, file = "NEUS40gfModelPerformance_v01.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(NEUS40gfa)[importance(NEUS40gfa) > 0])

#
# Overall Importance
png(file = "NEUS40Importance_v01.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(NEUS40gfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "NEUS40splitRatio_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(NEUS40gfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "NEUS40splitImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(NEUS40gfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "NEUS40splitData_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(NEUS40gfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="NEUS40cumImportanceSplit_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS40gfa, 
          imp.vars= imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 6,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = "NEUS40cumImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS40gfa, 
          imp.vars = imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 6,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(NEUS40gfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1974:2013)
NEUS40PCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="NEUS40PCA_v01.png", PCbiplot(NEUS40PCs),
       width = 83, height = 83, units = "mm", scale = 2)
###########

#find names of important variables
imp.vars
#pull data to get threshold quantiles from split density
a<-NEUS40gfa$d$Exploitation$x
b<-NEUS40gfa$d$Population$x
c<-NEUS40gfa$d$Landings$x
d<-NEUS40gfa$d$AMO$x
e<-NEUS40gfa$d$SST$x
f<-NEUS40gfa$d$Seafood$x
g<-NEUS40gfa$d$NAO_w$x
h<-NEUS40gfa$d$GS$x
i<-NEUS40gfa$d$MEI$x
NEUS40.dens<-cbind(a,b,c,d,e,f,g,h,i)

NEUS40.dens<-as.data.frame(NEUS40.dens)
colnames(NEUS40.dens)<-c("Exploitation", "Population", "Landings", "AMO", "SST", "Seafood", "NAO_w", "GS", "MEI")
NEUS40.d<-sapply(NEUS40.dens, summary)
write.csv(NEUS40.d, file="NEUS40.dens.split.csv")

########################
####50 year##########
######################


lev50 <- floor(log2(nrow(dat50) * 0.368/2))
#colnames(dat)[!colnames(dat) %in% c(dri.name, ind.name)]\
#dat<-data.frame(dat)
#
## GF analysis ##
NEUS50gfa <-   gradientForest(data = dat50, 
                              predictor.vars = dri.name50, 
                              response.vars = ind.name50,
                              ntree = 1000, 
                              transform = NULL,
                              maxLevel = lev50,
                              corr.threshold = 0.5, 
                              compact = F,
                              trace = T)

#save Rdata

save(NEUS50gfa, file = "NEUS50gfa_v01.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(NEUS50gfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- NEUS50gfa$imp.rsq
indicatorNA <- colnames(indicatorIMP)[apply(indicatorIMP, 2, function(x)all(is.na(x)))]
indicatorIMP <- indicatorIMP[apply(indicatorIMP, 2, function(x)!all(is.na(x))),]
#
indicatorDF <- data.frame("VARIABLE" = colnames(indicatorIMP),
                          "TYPE" = "INDICATOR", melt(apply(indicatorIMP, 2, mean, na.rm = T),
                                                     variable.name = "VARIABLE",
                                                     value.name = "MEAN"),
                          melt(apply(indicatorIMP, 2, min, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MIN"),
                          melt(apply(indicatorIMP, 2, max, na.rm = T),
                               variable.name = "VARIABLE", 
                               value.name = "MAX"))
#
pressureIMP <- NEUS50gfa$imp.rsq
pressureNA <- row.names(pressureIMP)[apply(pressureIMP, 1, function(x)all(is.na(x)))]
pressureIMP <- pressureIMP[apply(pressureIMP, 1, function(x)!all(is.na(x))),]
pressureDF <- data.frame("VARIABLE" = row.names(pressureIMP),
                         "TYPE" = "PRESSURE", melt(apply(pressureIMP, 1, mean, na.rm = T),
                                                   variable.name = "VARIABLE",
                                                   value.name = "MEAN"),
                         melt(apply(pressureIMP, 1, min, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MIN"),
                         melt(apply(pressureIMP, 1, max, na.rm = T),
                              variable.name = "VARIABLE", 
                              value.name = "MAX"))
#
indicatorDF[,c(3:5)] <- apply(indicatorDF[,c(3:5)], 2, signif, 2)
pressureDF[,c(3:5)] <- apply(pressureDF[,c(3:5)], 2, signif, 2)
# indicatorDF$VARIABLE <- factor(indicatorDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))
# pressureDF$VARIABLE <- factor(pressureDF$VARIABLE, labels = c(ENTER APPROPRIATE NAMES HERE))

#pressureDF <- rbind(pressureDF, cbind(VARIABLE = pressureNA,TYPE = "PRESSURE", MEAN = "NA", MIN = "NA", MAX = "NA"))
#
modPerformance <- rbind(pressureDF, indicatorDF)
##
write.csv(modPerformance, file = "NEUS50gfModelPerformance_v01.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(NEUS50gfa)[importance(NEUS50gfa) > 0])

#
# Overall Importance
png(file = "NEUS50Importance_v01.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(NEUS50gfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "NEUS50splitRatio_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(NEUS50gfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "NEUS50splitImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(NEUS50gfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "NEUS50splitData_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(NEUS50gfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="NEUS50cumImportanceSplit_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS50gfa, 
          imp.vars= imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 6,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = "NEUS50cumImportance_v01.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v01.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUS50gfa, 
          imp.vars = imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 6,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(NEUS50gfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1964:2013)
NEUS50PCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="NEUS50PCA_v01.png", PCbiplot(NEUS50PCs),
       width = 83, height = 83, units = "mm", scale = 2)
###########

#find names of important variables
imp.vars
#pull data to get threshold quantiles from split density
a<-NEUS50gfa$d$Exploitation$x
b<-NEUS50gfa$d$Population$x
c<-NEUS50gfa$d$Landings$x
d<-NEUS50gfa$d$AMO$x
e<-NEUS50gfa$d$SST$x
f<-NEUS50gfa$d$Seafood$x
g<-NEUS50gfa$d$NAO_w$x
h<-NEUS50gfa$d$GS$x
i<-NEUS50gfa$d$MEI$x
NEUS50.dens<-cbind(a,b,c,d,e,f,g,h,i)

NEUS50.dens<-as.data.frame(NEUS50.dens)
colnames(NEUS50.dens)<-c("Exploitation", "Population", "Landings", "AMO", "SST", "Seafood", "NAO_w", "GS", "MEI")
NEUS50.d<-sapply(NEUS50.dens, summary)
write.csv(NEUS50.d, file="NEUS50.dens.split.csv")