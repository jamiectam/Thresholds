set.seed(627)
#########
require(reshape2)
require(gradientForest)
require(ggplot2)
require(extendedForest)

###########
GOMenv1<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc")]
GOMenv2<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc")]
GOMenv3<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP")]
GOMenv4<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP")]
GOMenv5<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "Chlorophyll")]
GOMenv6<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "Chlorophyll")]
GOMenv7<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc", "GDP")]
GOMenv8<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP", "GDP.inc")]
GOMenv9<-GOMenv[,c("Landings", "Landings..1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater")]
GOMenv10<-GOMenv[,c("Landings", "Landings..1", "Exploitation", "Exploitation..1",
                      "SST", "AMO", "MEI", "Freshwater", "GDP")]
GOMenv11<-GOMenv[,c("Landings", "Exploitation",
                      "SST", "AMO", "MEI", "Freshwater", "GDP")]
GOMenv12<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom")]
GOMenv13<-GOMenv[,c("Landings", "Landings..1", "Exploitation", "Exploitation..1",
                      "SST", "MEI", "Freshwater.anom")]
GOMenv14<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings..1", "Exploitation", "Exploitation..1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
GOMenv15<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Exploitation",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
GOMenv12b<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                    "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom")]
GOMenv14a<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                    "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
GOMenv15<-GOMenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                    "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom", "Chlorophyll", "Currents", "AWP", "Hypoxic.area")]

#Set up data (use region specific env and bio dataframes)
dat<-cbind(GOMenv15, GOMbio)

ind.name<-colnames(GOMbio)
dri.name<-colnames(GOMenv15)

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
GOMEXgfa <-   gradientForest(data = dat, 
                             predictor.vars = dri.name, 
                             response.vars = ind.name,
                             ntree = 2000, 
                             transform = NULL,
                             maxLevel = lev,
                             corr.threshold = 0.5, 
                             compact = F,
                             trace = T)

#save Rdata

save(GOMEXgfa, file = "GOMEXgfa_v2.0_15.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(GOMEXgfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- GOMEXgfa$imp.rsq
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
pressureIMP <- GOMEXgfa$imp.rsq
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
write.csv(modPerformance, file = "GOMEXgfModelPerformance_v2.0_15.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(GOMEXgfa)[importance(GOMEXgfa) > 0])
#
# Overall Importance
png(file = "GOMEXImportance_v2.0_15.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(GOMEXgfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "GOMEXsplitRatio_v2.0_15.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "GOMEXE-splitRatio_v2.0_15.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(GOMEXgfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "GOMEXsplitImportance_v2.0_15.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "GOMEXE-splitImportance_v2.0_15.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(GOMEXgfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "GOMEXsplitData_v2.0_15.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "GOMEXE-splitData_v2.0_15.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(GOMEXgfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="GOMEXcumImportanceSplit_v2.0_15.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "GOMEXE-cumImportanceSplit_v2.0_15.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(GOMEXgfa, 
          imp.vars= imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 7,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Cumulative Importance (total)
png(file = "GOMEXcumImportance_v2.0_15.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "GOMEXE-cumImportance_v2.0_15.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(GOMEXgfa, 
          imp.vars = imp.vars, 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = FALSE,
          show.overall = TRUE,
          common.scale = TRUE,
          leg.nspecies = 7,
          cex.lab = 0.9, cex.legend = 0.7, cex.axis = 0.8, line.ylab = 0.8)
dev.off()
#
# Collect PCA info
Trns_grid <- predict(GOMEXgfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1992:2010)
GOMEXPCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="GOMEXPCA_v2.0_15.png", PCbiplot(GOMEXPCs),
       width = 83, height = 83, units = "mm", scale = 2)

########################

load("GOMEXgfa_v2.0_15.RDATA")
imp.vars

Population<-GOMEXgfa$d$Population.inc$x
SST<-GOMEXgfa$d$SST$x
AMO<-GOMEXgfa$d$AMO$x
Seafood<-GOMEXgfa$d$Seafood$x
Exploitation<-GOMEXgfa$d$Exploitation$x
Exploitation_1<-GOMEXgfa$d$Exploitation_1$x
Currents<-GOMEXgfa$d$Currents$x
Hypoxic.area<-GOMEXgfa$d$Hypoxic.area$x
AWP<-GOMEXgfa$d$AWP$x
Landings<-GOMEXgfa$d$Landings$x
Landings_1<-GOMEXgfa$d$Landings_1$x
MEI<-GOMEXgfa$d$MEI$x
Chlorophyll<-GOMEXgfa$d$Chlorophyll$x
GDP<-GOMEXgfa$d$GDP.inc$x
Freshwater<-GOMEXgfa$d$Freshwater.anom$x

GOMEX.dens<-cbind(Population, SST, AMO, Seafood, Exploitation, Currents,  Hypoxic.area, AWP, Landings, 
                   MEI, Landings_1, Exploitation_1, Chlorophyll, GDP, Freshwater)

GOMEX.dens<-as.data.frame(GOMEX.dens)

GOMEX.d<-sapply(GOMEX.dens, summary)
write.csv(GOMEX.d, file="GOMEX.dens.split.csv")