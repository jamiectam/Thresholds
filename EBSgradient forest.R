set.seed(627)
#########
require(reshape2)
require(gradientForest)
require(ggplot2)
require(extendedForest)

###########
EBSenv1<-EBSenv[,-14]

#Set up data (use region specific env and bio dataframes)
dat<-cbind(EBSenv1, EBSbio)

ind.name<-colnames(EBSbio)
dri.name<-colnames(EBSenv1)

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
EBSgfa <-   gradientForest(data = dat, 
                          predictor.vars = dri.name, 
                          response.vars = ind.name,
                          ntree = 1500, 
                          transform = NULL,
                          maxLevel = lev,
                          corr.threshold = 0.9, 
                          compact = F,
                          trace = T)

#save Rdata

save(EBSgfa, file = "EBSgfa_v2016_03.RDATA")



#
##########################
# Model Importance Table #
##########################
#
# Cumulative importance
var.order <- names(importance(EBSgfa, type = "Weighted", sort = TRUE))
#
indicatorIMP <- EBSgfa$imp.rsq
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
pressureIMP <- EBSgfa$imp.rsq
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
write.csv(modPerformance, file = "EBSgfModelPerformance_v2016_03.csv")
#         
#########
# PLOTS #
#########
#
imp.vars <- names(importance(EBSgfa)[importance(EBSgfa) > 0])
#
# Overall Importance
png(file = "EBSImportance_v2016_03.png", width = 83, height = 83, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.75, 6.5, 0.1, .5),
    omi = c(0, 0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
impPlot(EBSgfa,
        cex.main = 0.7)
mtext(expression(paste(R^2, " weighted importance")), side = 1, line = 1.75, cex = .75)
dev.off()
#
#
# Split Ratio
png(file = "EBSsplitRatio_v2016_03.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "EBSE-splitRatio_v2016_03.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spRatio(EBSgfa,
        imp.vars = imp.vars,
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.4, cex.axis = 0.6,
        cex.lab = 0.7, line.ylab = 0.9)
dev.off()

# Density of Splits
png(file = "EBSsplitImportance_v2016_03.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "EBSE-splitImportance_v2016_03.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spImportance(EBSgfa,
             imp.vars = imp.vars,
             #imp.vars.names = c("Coastal_engineering"),
             leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
             cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Density of Data
png(file = "EBSsplitData_v2016_03.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "EBSE-splitData_v2016_03.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(.75, 0.05, 0), 
    mar = c(2.25, 1, 0.1, .5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", 
    tck = -0.015)
spData(EBSgfa,
       imp.vars = imp.vars,
       #           imp.vars.names = #c(ADD NAMES HERE),
       leg.posn = "topright", cex.legend = 0.4, cex.axis = 0.6,
       cex.lab = 0.7, line.ylab = 0.9)
dev.off()
#
# Cumulative Importance (split)
png(file ="EBScumImportanceSplit_v2016_03.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "EBSE-cumImportanceSplit_v2016_03.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(EBSgfa, 
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
png(file = "EBScumImportance_v2016_03.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "EBSE-cumImportance_v2016_03.tiff"), width = 173.5, height = 173.5, units = "mm", res = 600)
#
par(mgp = c(1.0, 0.1, 0),
    mar = c(2.25, 1, 0.1, 0.5),
    omi = c(0,0.3, 0.1, 0),
    family = "sans", tck = -0.015)
spCumPlot(EBSgfa, 
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
Trns_grid <- predict(EBSgfa, dat[, imp.vars])
row.names(Trns_grid) <- c(1982:2013)
EBSPCs <- prcomp(Trns_grid[, imp.vars])

#add this to function
require(grid)

ggsave(file="EBSPCA_v2016_03.png", PCbiplot(EBSPCs),
       width = 83, height = 83, units = "mm", scale = 2)

###############
imp.vars

Population<-EBSgfa$d$Population.inc$x
SST<-EBSgfa$d$SST$x
PDO<-EBSgfa$d$PDO$x
Seafood<-EBSgfa$d$Seafood$x
Exploitation<-EBSgfa$d$Exploitation$x
Max.Trawled.Area<-EBSgfa$d$Max.Trawled.Area$x
Cold.pool<-EBSgfa$d$Cold.pool$x
Ice.Retreat<-EBSgfa$d$Ice.Retreat$x
BTMP<-EBSgfa$d$BTMP$x
Landings<-EBSgfa$d$Landings$x
MEI<-EBSgfa$d$MEI$x
NOD<-EBSgfa$d$NOD$x
NPI<-EBSgfa$d$NPI$x

EBS.dens<-cbind(MEI, PDO, NOD, Exploitation, Landings, Population, 
                Seafood, Ice.Retreat, NPI, Max.Trawled.Area, Cold.pool,
                BTMP, SST) 

EBS.dens<-as.data.frame(EBS.dens)

EBS.d<-sapply(EBS.dens, summary)
write.csv(EBS.d, file="EBS.dens.split.csv")