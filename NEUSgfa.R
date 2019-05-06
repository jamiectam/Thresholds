set.seed(88)
#########
require(reshape2)
require(gradientForest)
require(ggplot2)
require(extendedForest)
require(grid)

NEUSbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/NEUSbio.csv")
NEUSenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/NEUSenv.csv", header=T)
EBSbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/EBSbio.csv")
EBSenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/EBSenv.csv")
ccbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/ccbio.csv")
ccenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/ccenv.csv")
GOMbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/GOMbio.csv")
GOMenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/GOMenv.csv")
PIbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/PIbio.csv")
PIenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/PIenv.csv")
###########

NEUSenv1<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings.1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc")]
NEUSenv2<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc")]
NEUSenv3<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings.1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP")]
NEUSenv4<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP")]
NEUSenv5<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings.1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "Chlorophyll")]
NEUSenv6<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "Chlorophyll")]
NEUSenv7<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings.1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc", "GDP")]
NEUSenv8<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP", "GDP.inc")]
NEUSenv9<-NEUSenv[,c("Landings", "Landings.1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater")]
NEUSenv10<-NEUSenv[,c("Landings", "Landings.1", "Exploitation", "Exploitation..1",
                     "SST", "AMO", "MEI", "Freshwater", "GDP")]
NEUSenv11<-NEUSenv[,c("Landings", "Exploitation",
                      "SST", "AMO", "MEI", "Freshwater", "GDP")]
NEUSenv12<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings.1", "Exploitation", "Exploitation..1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom")]
NEUSenv13<-NEUSenv[,c("Landings", "Landings.1", "Exploitation", "Exploitation..1",
                      "SST", "MEI", "Freshwater.anom")]
NEUSenv14<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings.1", "Exploitation", "Exploitation..1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
NEUSenv15<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
NEUSenv00sl<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation", 
                     "SST", "AMO", "MEI", "Freshwater", "GDP.inc")]
NEUSenv12b<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom")]
NEUSenv12c<-NEUSenv[,c("Population.inc", "Seafood", "Landings0", "Exploitation",
                       "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom")]
NEUSenv12d<-NEUSenv[,c("Population.inc", "Seafood", "Landings0", "Landings.10", "Exploitation", "Exploitation_1",
                       "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom")]
NEUSenv14a<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]

NEUSenv14b<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Exploitation",
                       "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
NEUSenv14c<-NEUSenv[,c("Population.inc", "Seafood", "Landings0", "Landings.10", "Exploitation", "Exploitation_1",
                       "SST", "AMO", "MEI", "GDP.inc", "Freshwater.total.anom")]
NEUSenv15<-NEUSenv[,c("Population.inc", "Seafood", "Landings", "Landings_1", "Exploitation", "Exploitation_1",
                      "SST", "AMO", "MEI", "GDP.inc", "Freshwater.anom", "NAO_w", "Wind", "Chlorophyll", "GS")]
#Set up data (use region specific env and bio dataframes)
dat<-cbind(NEUSenv15, NEUSbio)

ind.name<-colnames(NEUSbio)
dri.name<-colnames(NEUSenv15)

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
NEUSgfa <-   gradientForest(data = dat, 
                            predictor.vars = dri.name, 
                            response.vars = ind.name,
                            ntree = 2000, 
                            transform = NULL,
                            maxLevel = lev,
                            corr.threshold = 0.98, 
                            compact = F,
                            trace = T)

#save Rdata

save(NEUSgfa, file = "NEUSgfa_v2.0_15a.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(NEUSgfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- NEUSgfa$imp.rsq
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
pressureIMP <- NEUSgfa$imp.rsq
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
write.csv(modPerformance, file = "NEUSgfModelPerformance_v2.0_15a.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(NEUSgfa)[importance(NEUSgfa) > 0])

#
# Overall Importance
png(file = "NEUSImportance_v2.0_15a.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(NEUSgfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "NEUSsplitRatio_v2.0_15a.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitRatio_v2.0_15a.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(NEUSgfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "NEUSsplitImportance_v2.0_15a.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitImportance_v2.0_15a.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(NEUSgfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "NEUSsplitData_v2.0_15a.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-splitData_v2.0_15a.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(NEUSgfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="NEUScumImportanceSplit_v2.0_15a.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportanceSplit_v2.0_15a.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUSgfa, 
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
png(file = "NEUScumImportance_v2.0_15a.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "CCE-cumImportance_v2.0_15a.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(NEUSgfa, 
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
Trns_grid <- predict(NEUSgfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1964:2013)
NEUSPCs <- prcomp(Trns_grid[, imp.vars])
#
ggsave(file="NEUSPCA_v2.0_15a.png", PCbiplot(NEUSPCs),
       width = 83, height = 83, units = "mm", scale = 2)
###########

load("NEUSgfa_v2.0_15a.RDATA")
#find names of important variables
#imp.vars
#pull data to get threshold quantiles from split density

Population<-NEUSgfa$d$Population.inc$x
SST<-NEUSgfa$d$SST$x
Seafood<-NEUSgfa$d$Seafood$x
Exploitation<-NEUSgfa$d$Exploitation$x
Exploitation_1<-NEUSgfa$d$Exploitation_1$x
AMO<-NEUSgfa$d$AMO$x
MEI<-NEUSgfa$d$MEI$x
Landings<-NEUSgfa$d$Landings$x
Landings_1<-NEUSgfa$d$Landings_1$x
Wind<-NEUSgfa$d$Wind$x
NAO_w<-NEUSgfa$d$NAO_w$x
GS<-NEUSgfa$d$GS$x
GDP<-NEUSgfa$d$GDP.inc$x
Freshwater<-NEUSgfa$d$Freshwater.anom$x
Chlorophyll<-NEUSgfa$d$Chlorophyll$x

NEUS.dens<-cbind(Population, SST, AMO, Seafood, Exploitation, MEI,  Landings, 
                 NAO_w, GS, GDP, Freshwater, Landings_1, Exploitation_1, Wind, Chlorophyll)

NEUS.dens<-as.data.frame(NEUS.dens)

neus.d<-sapply(NEUS.dens, summary)
write.csv(neus.d, file="NEUS.dens.split1.csv")
