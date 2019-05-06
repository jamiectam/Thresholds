set.seed(627)
#########
require(reshape2)
require(gradientForest)
require(ggplot2)
require(extendedForest)
require(grid)

###########
ccenv1<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "PDO", "MEI", "Freshwater", "GDP.inc")]
ccenv2<-ccenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "PDO", "MEI", "Freshwater", "GDP.inc")]
ccenv3<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "PDO", "MEI", "Freshwater", "GDP")]
ccenv4<-ccenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "PDO", "MEI", "Freshwater", "GDP")]
ccenv5<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "PDO", "MEI", "Freshwater", "Chl..µg.l.")]
ccenv6<-ccenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "PDO", "MEI", "Freshwater", "Chl..µg.l.")]
ccenv7<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                   "SST", "PDO", "MEI", "Freshwater", "GDP.inc", "GDP")]
ccenv8<-ccenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                   "SST", "PDO", "MEI", "Freshwater", "GDP", "GDP.inc")]
ccenv9<-ccenv[,c("Landings", "Landings..1", "Exploitation", "Exploitation..1",
                   "SST", "PDO", "MEI", "Freshwater")]
ccenv10<-ccenv[,c("Landings", "Landings..1", "Exploitation", "Exploitation..1",
                    "SST", "PDO", "MEI", "Freshwater", "GDP")]
ccenv11<-ccenv[,c("Landings", "Exploitation",
                    "SST", "PDO", "MEI", "Freshwater", "GDP")]
ccenv12<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                    "SST", "PDO", "MEI", "GDP.inc", "Freshwater.anom")]
ccenv13<-ccenv[,c("Landings", "Landings_1", "Exploitation", "Exploitation_1",
                    "SST", "MEI", "Freshwater.anom")]
ccenv14<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                    "SST", "PDO", "MEI", "GDP.inc", "Freshwater.total.anom")]
ccenv15<-ccenv[,c("Population.inc", "Seafood", "Landings", "Exploitation",
                    "SST", "PDO", "MEI", "GDP.inc", "Freshwater.total.anom")]
ccenv12b<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                     "SST", "PDO", "MEI", "GDP.inc", "Freshwater.anom")]
ccenv14a<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                     "SST", "PDO", "MEI", "GDP.inc", "Freshwater.total.anom")]
ccenv15<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                  "SST", "PDO", "MEI", "GDP.inc", "NPGO", "Chlorophyll")]
ccenv15b<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                  "SST", "PDO", "MEI", "GDP.inc", "NPGO", "Chlorophyll", "Freshwater.anom")]
ccenv15c<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                   "SST", "PDO", "MEI", "GDP.inc", "NPGO", "Chlorophyll", "Freshwater.anom", "TUMI")]
ccenv15d<-ccenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                   "SST", "PDO", "MEI", "GDP.inc", "NPGO", "Chlorophyll", "Freshwater.anom", "TUMI", "Sea.Level")]



#Set up data (use region specific env and bio dataframes)
dat<-cbind(ccenv15d, ccbio)

ind.name<-colnames(ccbio)
dri.name<-colnames(ccenv15d)

# Many of the pressure variables do not have full data (some NA):
## Option 1- Impute:
dat <- na.roughfix(dat)
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
lev <- floor(log2(nrow(dat) * 0.368/2))
#colnames(dat)[!colnames(dat) %in% c(dri.name, ind.name)]\
#dat<-data.frame(dat)
#
## GF analysis ##
ccgfa <-   gradientForest(data = dat, 
                          predictor.vars = dri.name, 
                          response.vars = ind.name,
                          ntree = 1500, 
                          transform = NULL,
                          maxLevel = lev,
                          corr.threshold = 0.5, 
                          compact = F,
                          trace = T)

#save Rdata

save(ccgfa, file = "ccgfa_v2.0_15d.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(ccgfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- ccgfa$imp.rsq
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
pressureIMP <- ccgfa$imp.rsq
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
write.csv(modPerformance, file = "ccgfModelPerformance_v2.0_15d.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(ccgfa)[importance(ccgfa) > 0])
#
# Overall Importance
png(file = "ccImportance_v2.0_15d.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(ccgfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "ccsplitRatio_v2.0_15d.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v2.0_15d.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(ccgfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "ccsplitImportance_v2.0_15d.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v2.0_15d.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(ccgfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "ccsplitData_v2.0_15d.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v2.0_15d.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(ccgfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="cccumImportanceSplit_v2.0_15d.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v2.0_15d.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)

spCumPlot(ccgfa, 
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
png(file = "cccumImportance_v2.0_15d.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v2.0_15d.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(ccgfa, 
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
Trns_grid <- predict(ccgfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1981:2012)
ccPCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="ccPCA_v2.0_15d.png", PCbiplot(ccPCs),
       width = 83, height = 83, units = "mm", scale = 2)

##################

load("ccgfa_v2.0_15c.RDATA")
imp.vars

Population<-ccgfa$d$Population.inc$x
SST<-ccgfa$d$SST$x
PDO<-ccgfa$d$PDO$x
Seafood<-ccgfa$d$Seafood$x
Exploitation<-ccgfa$d$Exploitation$x
Exploitation_1<-ccgfa$d$Exploitation_1$x
MEI<-ccgfa$d$MEI$x
GDP<-ccgfa$d$GDP.inc$x
Chlorophyll<-ccgfa$d$Chlorophyll$x
Landings<-ccgfa$d$Landings$x
Landings_1<-ccgfa$d$Landings_1$x
Freshwater<-ccgfa$d$Freshwater.anom$x
NPGO<-ccgfa$d$NPGO$x
TUMI<-ccgfa$d$TUMI$x

#Commercial.shipping<-ccgfa$d$Commercial.shipping$x
#Atmospheric.pollution<-ccgfa$d$Atmospheric.pollution$x
#Habitat_modification<-ccgfa$d$Habitat_modification$x

cc.dens<-cbind(Population, SST, PDO, Seafood, Exploitation, MEI, TUMI, Landings, Exploitation_1, Landings_1, GDP, Freshwater, 
               Chlorophyll, NPGO)
  
cc.dens<-as.data.frame(cc.dens)

cc.d<-sapply(cc.dens, summary)
write.csv(cc.d, file="cc.dens.split.csv")