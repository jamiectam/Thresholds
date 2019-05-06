require(ggplot2)
require(grid)
require(gridExtra)
require(png)

load("NEUSgfa_v2016_04.RDATA")
load("ccgfa_v2016_04.RDATA")
load("EBSgfa_v2016_03.RDATA")
load("GOMEXgfa_v2016_2005.RDATA")
load("PIgfa_v04.RDATA")


###############importance plots###################

EBSimp<-as.raster(readPNG("EBSImportance_v2016_03.png"))
ccimp<-as.raster(readPNG("ccImportance_v2016_04.png"))
NEUSimp<-as.raster(readPNG("NEUSImportance_v2016_04.png"))
GOMEXimp<-as.raster(readPNG("GOMEXImportance_v2016_2005.png"))
WHimp<-as.raster(readPNG("PIImportance_v04.png"))

EBSimp<-rasterGrob(EBSimp, interpolate=F)
ccimp<-rasterGrob(ccimp, interpolate=F)
NEUSimp<-rasterGrob(NEUSimp, interpolate=F)
GOMEXimp<-rasterGrob(GOMEXimp, interpolate=F)
WHimp<-rasterGrob(WHimp, interpolate=F)

png(file ="importanceplots_v04.png", width = 600, height = 800, units = "mm", res = 600)
grid.arrange(EBSimp, ccimp, NEUSimp, GOMEXimp, WHimp, nrow=3, ncol=2)
grid.text("a)", x=unit(0.10, "npc"), y=unit(0.97, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("b)", x=unit(0.65, "npc"), y=unit(0.97, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("c)", x=unit(0.10, "npc"), y=unit(0.65, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("d)", x=unit(0.65, "npc"), y=unit(0.65, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("e)", x=unit(0.10, "npc"), y=unit(0.3, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))

dev.off()



#######Spcumplot##############
png(file ="SpCompPDOAMO_v05.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mfrow=c(2,2))

   # mgp = c(1.0, 0.1, 0),
   # mar = c(2.25, 1, 0.1, 0.5),
    #omi = c(0,0.3, 0.1, 0),
   # family = "sans", tck = -0.015)

spCompPlot(EBSgfa, 
          imp.vars= c("PDOa"), 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 7,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)

spCompPlot(ccgfa, 
          imp.vars= c("PDO"), 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 7,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)

spCompPlot(NEUSgfa, 
          imp.vars= c("AMO"), 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 7,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)

spCompPlot(GOMEXgfa, 
          imp.vars= c("AMO"), 
          #           imp.vars.names = #c(ADD NAMES HERE),
          show.species = TRUE,
          legend = TRUE,
          show.overall = FALSE,
          common.scale = TRUE,
          leg.nspecies = 7,
          #           leg.posn = "topright",
          cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()

#######################
#cumulation plots total
######################

png(file ="totalCompPDOAMO_v05.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mfrow=c(2,2))

# mgp = c(1.0, 0.1, 0),
# mar = c(2.25, 1, 0.1, 0.5),
#omi = c(0,0.3, 0.1, 0),
# family = "sans", tck = -0.015)

spCompPlot(EBSgfa, 
           imp.vars= c("PDOa"), 
           #           imp.vars.names = #c(ADD NAMES HERE),
           show.species = FALSE,
           legend = TRUE,
           show.overall = TRUE,
           common.scale = TRUE,
           leg.nspecies = 7,
           #           leg.posn = "topright",
           cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)

spCompPlot(ccgfa, 
           imp.vars= c("PDO"), 
           #           imp.vars.names = #c(ADD NAMES HERE),
           show.species = FALSE,
           legend = TRUE,
           show.overall = TRUE,
           common.scale = TRUE,
           leg.nspecies = 7,
           #           leg.posn = "topright",
           cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)

spCompPlot(NEUSgfa, 
           imp.vars= c("AMO"), 
           #           imp.vars.names = #c(ADD NAMES HERE),
           show.species = FALSE,
           legend = TRUE,
           show.overall = TRUE,
           common.scale = TRUE,
           leg.nspecies = 7,
           #           leg.posn = "topright",
           cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)

spCompPlot(GOMEXgfa, 
           imp.vars= c("AMO"), 
           #           imp.vars.names = #c(ADD NAMES HERE),
           show.species = FALSE,
           legend = TRUE,
           show.overall = TRUE,
           common.scale = TRUE,
           leg.nspecies = 7,
           #           leg.posn = "topright",
           cex.lab = 0.9, cex.legend = 1, cex.axis = 0.8, line.ylab = 0.8)
dev.off()

########comparative thresholds######
EBS.d<-read.csv("EBS.dens.split.csv")
cc.d<-read.csv("cc.dens.split.csv")
GOMEX.d<-read.csv("GOMEX.dens.split.csv")
NEUS.d<-read.csv("NEUS.dens.split.csv")
PI.d<-read.csv("PI.dens.split.csv")


EBS.a<-EBS.d[c(2, 4, 5),c("Population", "Seafood", "Landings", "Exploitation", "SST", "MEI", "PDO")]
cc.a<-cc.d[c(2, 4, 5),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "PDO")]
GOMEX.a<-GOMEX.d[c(2, 4, 5),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "AMO")]
NEUS.a<-NEUS.d[c(2, 4, 5),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "AMO")]
PI.a<-PI.d[c(2, 4, 5),c("Population", "Seafood", "Landings", "SST",  "MEI", "PDO")]

IEA<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
IEA2<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US")

Population<-cbind(EBS.a$Population, cc.a$Population, GOMEX.a$Population, NEUS.a$Population, PI.a$Population)
colnames(Population)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Population)<-c("Upper", "Mean", "Lower")
Population<-t(Population)
Population<-as.data.frame(Population)
Population<-cbind(Population, IEA)


SST<-cbind(EBS.a$SST, cc.a$SST, GOMEX.a$SST, NEUS.a$SST, PI.a$SST)
colnames(SST)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(SST)<-c("Upper", "Mean", "Lower")
SST<-t(SST)
SST<-as.data.frame(SST)
SST<-cbind(SST, IEA)

MEI<-cbind(EBS.a$MEI, cc.a$MEI, GOMEX.a$MEI, NEUS.a$MEI, PI.a$MEI)
colnames(MEI)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(MEI)<-c("Upper", "Mean", "Lower")
MEI<-t(MEI)
MEI<-as.data.frame(MEI)
MEI<-cbind(MEI, IEA)

Seafood<-cbind(EBS.a$Seafood, cc.a$Seafood, GOMEX.a$Seafood, NEUS.a$Seafood, PI.a$Seafood)
colnames(Seafood)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Seafood)<-c("Upper", "Mean", "Lower")
Seafood<-t(Seafood)
Seafood<-as.data.frame(Seafood)
Seafood<-cbind(Seafood, IEA)

AMOPDO<-cbind(EBS.a$PDO, cc.a$PDO, GOMEX.a$AMO, NEUS.a$AMO, PI.a$PDO)
colnames(AMOPDO)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(AMOPDO)<-c("Upper", "Mean", "Lower")
AMOPDO<-t(AMOPDO)
AMOPDO<-as.data.frame(AMOPDO)
AMOPDO<-cbind(AMOPDO, IEA)

Landings<-cbind(EBS.a$Landings, cc.a$Landings, GOMEX.a$Landings, NEUS.a$Landings, PI.a$Landings)
colnames(Landings)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Landings)<-c("Upper", "Mean", "Lower")
Landings<-t(Landings)
Landings<-as.data.frame(Landings)
Landings<-cbind(Landings, IEA)

Exploitation<-cbind(EBS.a$Exploitation, cc.a$Exploitation, GOMEX.a$Exploitation, NEUS.a$Exploitation)
colnames(Exploitation)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US")
rownames(Exploitation)<-c("Upper", "Mean", "Lower")
Exploitation<-t(Exploitation)
Exploitation<-as.data.frame(Exploitation)
Exploitation<-cbind(Exploitation, IEA2)

# EBS.b<-EBS.d[c(2, 3, 4, 5),c("SST","PDO")]
# cc.b<-cc.d[c(2, 3, 4, 5),c( "SST", "PDO")]
# GOMEX.b<-GOMEX.d[c(2, 3, 4, 5),c("SST", "AMO")]
# NEUS.b<-NEUS.d[c(2, 3, 4, 5),c("SST", "AMO")]
# 
# all.dens<-rbind(EBS.a, cc.a, GOMEX.a, NEUS.a)
# 
# write.csv(all.dens, file="all.dens.csv")

# SST<-all.dens[,c(1,2,7)]
# SST<-dcast(SST, IEA~Stat)
# 
# 
#IEA<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
# 
# Population<-read.csv("PopulationThres.csv")
# Seafood<-read.csv("SeafoodThres.csv")
# Landings<-read.csv("LandingsThres.csv")
# Exploitation<-read.csv("ExploitationThres.csv")
#this defines the elements to go in the plot, both the x and y and upper and lower CIs



a<-ggplot(Population,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
        legend.position="none" )+
        ylim(0,5)+
        #labs(y=expression(paste("Population", " ", "(Individuals"," ", km^{-2},")", sep="")))
        ylab("Population increase (rate of change)")

b<-ggplot(Seafood,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
        legend.position="none" )+
        ylab("Seafood (lbs/capita USA)")

c<-ggplot(Landings,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
        legend.position="none" )+
        expand_limits(y=c(2,7))+
  labs(y=expression(paste("Landings", " ", "(t"," ", km^{-2},")", sep="")))
         

d<-ggplot(Exploitation,aes(x=IEA2,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
               legend.position="none" )+
        ylim(0,1.30)+
         ylab("Exploitation (Total landings/Total biomass)")

##############
# SST<-read.csv("SSTthresh.csv")
# PDOAMO<-read.csv("PDOAMOThres1.csv")

e<-ggplot(SST,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
        legend.position="none" )+
  ylab("SST (degree C)")

f<-ggplot(AMOPDO,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
        legend.position="none" )+
  ylab("AMO and PDO")


g<-ggplot(MEI,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
  geom_pointrange(size=1, colour='black')+
  coord_flip()+
  theme_bw()+
  theme(axis.line=element_line(colour='black'), 
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=16, colour='black'),
        axis.text.y=element_text(size=16, colour='black'),
        axis.title.x=element_text(size=14, colour='black'),
        axis.title.y=element_blank(),
        legend.position="none" )+
  ylab("MEI")


-----------
  #combine plots
-------------------

png(file ="AllThresh_v05.png", width = 500, height = 300, units = "mm", res = 600)
# # tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)
# 
grid.arrange(a,b,c,d,e,f,g, nrow=3, ncol=3)

# 
dev.off()


png(file ="AnthroThresh_v01.png", width = 500, height = 300, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

grid.arrange(a,b,c,d,nrow=2, ncol=2)

dev.off()  


png(file ="ClimateThresh_v01.png", width = 180, height = 300, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

grid.arrange(e,f,g, nrow=3, ncol=1)

dev.off()


---------------
  #single plots
  ---------------
################indv graphs

png(file ="Population_v03.png", width = 300, height = 150, units = "mm", res = 600)
a
dev.off()

png(file ="Seafood_v03.png", width = 300, height = 150, units = "mm", res = 600)
b
dev.off()

png(file ="Landings_v03.png", width = 300, height = 150, units = "mm", res = 600)
c
dev.off()

png(file ="Exploitation_v03.png", width = 300, height = 150, units = "mm", res = 600)
d
dev.off()

png(file ="SST_v03.png", width = 300, height = 150, units = "mm", res = 600)
e
dev.off()

png(file ="AMOPDO_v03.png", width = 300, height = 150, units = "mm", res = 600)
f
dev.off()

png(file ="MEI_v03.png", width = 300, height = 150, units = "mm", res = 600)
g
dev.off()

# Split Ratio

png(file ="CompSplitRatio.png", width = 173.5, height = 173.5, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

par(mfrow=c(2,2))

spRatio(EBSgfa,
        imp.vars = c("SST"),
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.9, cex.axis = 0.9)


spRatio(ccgfa,
        imp.vars = c("SST"),
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.9, cex.axis = 0.9,
        cex.lab = 0.9, line.ylab = 0.9)


spRatio(NEUSgfa,
        imp.vars = c("SST"),
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.9, cex.axis = 0.9,
        cex.lab = 0.9, line.ylab = 0.9)


spRatio(GOMEXgfa,
        imp.vars = c("SST"),
        #           imp.vars.names = #c(ADD NAMES HERE),
        leg.posn = "topright",
        cex.legend = 0.9, cex.axis = 0.9,
        cex.lab = 0.9, line.ylab = 0.9)
dev.off()
