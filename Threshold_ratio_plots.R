require(ggplot2)
require(reshape2)
require(grid)

EBS.d<-read.csv("EBS.dens.splitb.csv")
cc.d<-read.csv("cc.dens.splitb.csv")
GOMEX.d<-read.csv("GOMEX.dens.splitb.csv")
NEUS.d<-read.csv("NEUS.dens.split1b.csv")
PI.d<-read.csv("PI.dens.split1b.csv")


EBS.a<-EBS.d[c(2, 4, 5,7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST", "MEI", "PDO", "Exploitation_1", "Landings_1",
                                "Chlorophyll", "Freshwater", "GDP")]
cc.a<-cc.d[c(2, 4, 5,7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "PDO", "Exploitation_1", "Landings_1",
                              "Chlorophyll", "Freshwater", "GDP")]
GOMEX.a<-GOMEX.d[c(2, 4, 5,7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "AMO", "Exploitation_1", "Landings_1",
                                    "Chlorophyll", "Freshwater", "GDP")]
NEUS.a<-NEUS.d[c(2, 4, 5,7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "AMO", "Exploitation_1", "Landings_1",
                                  "Chlorophyll", "Freshwater", "GDP")]
PI.a<-PI.d[c(2, 4, 5,7,8,9),c("Population", "Seafood", "Landings", "SST",  "MEI", "PDO", "Landings_1", "GDP", "Chlorophyll")]



EBS.a.thresh<-EBS.d[c(7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST", "MEI", "PDO", "Exploitation_1", "Landings_1",
                               "Chlorophyll", "Freshwater", "GDP", "IEA")]
cc.a.thresh<-cc.d[c(7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "PDO", "Exploitation_1", "Landings_1",
                             "Chlorophyll", "Freshwater", "GDP", "IEA")]
GOMEX.a.thresh<-GOMEX.d[c(7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "AMO", "Exploitation_1", "Landings_1",
                                   "Chlorophyll", "Freshwater", "GDP", "IEA")]
NEUS.a.thresh<-NEUS.d[c(7,8,9),c("Population", "Seafood", "Landings", "Exploitation", "SST",  "MEI", "AMO", "Exploitation_1", "Landings_1",
                                 "Chlorophyll", "Freshwater", "GDP", "IEA")]
PI.a.thresh<-PI.d[c(7,8,9),c("Population", "Seafood", "Landings", "SST",  "MEI", "PDO", "Landings_1", "GDP", "Chlorophyll", "IEA")]


IEA<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawai'i")
IEA2<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US")

Population<-cbind(EBS.a.thresh$Population, cc.a.thresh$Population, GOMEX.a.thresh$Population, NEUS.a.thresh$Population, PI.a.thresh$Population)
colnames(Population)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Population)<-c( "Threshold", "Threshold", "Threshold")
Population<-t(Population)
Population<-as.data.frame(Population)
Population<-cbind(Population, IEA)


SST<-cbind(EBS.a$SST, cc.a$SST, GOMEX.a$SST, NEUS.a$SST, PI.a$SST)
colnames(SST)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(SST)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
SST<-t(SST)
SST<-as.data.frame(SST)
SST<-cbind(SST, IEA)

MEI<-cbind(EBS.a$MEI, cc.a$MEI, GOMEX.a$MEI, NEUS.a$MEI, PI.a$MEI)
colnames(MEI)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(MEI)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
MEI<-t(MEI)
MEI<-as.data.frame(MEI)
MEI<-cbind(MEI, IEA)

Seafood<-cbind(EBS.a$Seafood, cc.a$Seafood, GOMEX.a$Seafood, NEUS.a$Seafood, PI.a$Seafood)
colnames(Seafood)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Seafood)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Seafood<-t(Seafood)
Seafood<-as.data.frame(Seafood)
Seafood<-cbind(Seafood, IEA)

AMOPDO<-cbind(EBS.a$PDO, cc.a$PDO, GOMEX.a$AMO, NEUS.a$AMO, PI.a$PDO)
colnames(AMOPDO)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(AMOPDO)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
AMOPDO<-t(AMOPDO)
AMOPDO<-as.data.frame(AMOPDO)
AMOPDO<-cbind(AMOPDO, IEA)

Landings<-cbind(EBS.a$Landings, cc.a$Landings, GOMEX.a$Landings, NEUS.a$Landings, PI.a$Landings)
colnames(Landings)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Landings)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Landings<-t(Landings)
Landings<-as.data.frame(Landings)
Landings<-cbind(Landings, IEA)

Exploitation<-cbind(EBS.a$Exploitation, cc.a$Exploitation, GOMEX.a$Exploitation, NEUS.a$Exploitation)
colnames(Exploitation)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US")
rownames(Exploitation)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Exploitation<-t(Exploitation)
Exploitation<-as.data.frame(Exploitation)
Exploitation<-cbind(Exploitation, IEA2)

Landings_1<-cbind(EBS.a$Landings, cc.a$Landings_1, GOMEX.a$Landings_1, NEUS.a$Landings_1, PI.a$Landings_1)
colnames(Landings_1)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Landings_1)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Landings_1<-t(Landings_1)
Landings_1<-as.data.frame(Landings_1)
Landings_1<-cbind(Landings_1, IEA)

Exploitation_1<-cbind(EBS.a$Exploitation_1, cc.a$Exploitation_1, GOMEX.a$Exploitation_1, NEUS.a$Exploitation_1)
colnames(Exploitation_1)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US")
rownames(Exploitation_1)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Exploitation_1<-t(Exploitation_1)
Exploitation_1<-as.data.frame(Exploitation_1)
Exploitation_1<-cbind(Exploitation_1, IEA2)

Chlorophyll<-cbind(EBS.a$Chlorophyll, cc.a$Chlorophyll, GOMEX.a$Chlorophyll, NEUS.a$Chlorophyll, PI.a$Chlorophyll)
colnames(Chlorophyll)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(Chlorophyll)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Chlorophyll<-t(Chlorophyll)
Chlorophyll<-as.data.frame(Chlorophyll)
Chlorophyll<-cbind(Chlorophyll, IEA)

GDP<-cbind(EBS.a$GDP, cc.a$GDP, GOMEX.a$GDP, NEUS.a$GDP, PI.a$GDP)
colnames(GDP)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US", "West Hawaii")
rownames(GDP)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
GDP<-t(GDP)
GDP<-as.data.frame(GDP)
GDP<-cbind(GDP, IEA)

Freshwater<-cbind(EBS.a$Freshwater, cc.a$Freshwater, GOMEX.a$Freshwater, NEUS.a$Freshwater)
colnames(Freshwater)<-c("Bering Sea", "California Current", "Gulf of Mexico", "Northeast US")
rownames(Freshwater)<-c("Upper", "Mean", "Lower", "Threshold1", "Threshold2", "Threshold3")
Freshwater<-t(Freshwater)
Freshwater<-as.data.frame(Freshwater)
Freshwater<-cbind(Freshwater, IEA2)


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



a<-ggplot(Population,aes(x=IEA,y=Threshold,ymax=Upper,ymin=Lower))+
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
  ylim(0,1.0)+
  ylab("Exploitation (total landings/total biomass)")

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
  ylim(5,27)+
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

h<-ggplot(Landings_1,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
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
  expand_limits(y=c(2,6))+
  labs(y=expression(paste("Landings_1", " ", "(t"," ", km^{-2},")", sep="")))

i<-ggplot(Exploitation_1,aes(x=IEA2,y=Mean,ymax=Upper,ymin=Lower))+
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
  expand_limits(y=c(0,1.0))+
  ylab("Exploitation_1(total landings/total biomass)")

j<-ggplot(GDP,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
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
  ylab("GDP increase (USD)")

k<-ggplot(Freshwater,aes(x=IEA2,y=Mean,ymax=Upper,ymin=Lower))+
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
  ylab("Freshwater anomoalies")

l<-ggplot(Chlorophyll,aes(x=IEA,y=Mean,ymax=Upper,ymin=Lower))+
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
  labs(y=expression(paste("Chlorophyll", " ", "(mg"," ", m^{-3},")", sep="")))




-----------
  #combine plots
  -------------------
  
  png(file ="AllThresh2.0_v05.png", width = 500, height = 400, units = "mm", res = 600)
# # tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)
# 
grid.arrange(c,h,d,i,a,b,j,e,f,g,k,l, nrow=4, ncol=3)
grid.text("a)", x=unit(0.08, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("b)", x=unit(0.40, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("c)", x=unit(0.75, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("d)", x=unit(0.08, "npc"), y=unit(0.74, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("e)", x=unit(0.40, "npc"), y=unit(0.74, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("f)", x=unit(0.75, "npc"), y=unit(0.74, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("g)", x=unit(0.08, "npc"), y=unit(0.49, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("h)", x=unit(0.40, "npc"), y=unit(0.49, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("i)", x=unit(0.75, "npc"), y=unit(0.49, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("j)", x=unit(0.08, "npc"), y=unit(0.24, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("k)", x=unit(0.40, "npc"), y=unit(0.24, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("l)", x=unit(0.75, "npc"), y=unit(0.24, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))

# 
dev.off()


png(file ="AnthroThresh2.0_v04.png", width = 500, height = 300, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

grid.arrange(c,d,h,i,a,b,j,nrow=3, ncol=3)
grid.text("a)", x=unit(0.08, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("b)", x=unit(0.40, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("c)", x=unit(0.75, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("d)", x=unit(0.08, "npc"), y=unit(0.65, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("e)", x=unit(0.40, "npc"), y=unit(0.65, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("f)", x=unit(0.75, "npc"), y=unit(0.65, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("g)", x=unit(0.08, "npc"), y=unit(0.33, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
dev.off()  


png(file ="ClimateThresh2.0_v04.png", width = 500, height = 250, units = "mm", res = 600)
# tiff(file = paste0(figure.dir, "PIE-cumImportanceSplit_v05.tiff"), width = 173.5, height = 173.5, units = "mm", res = 1000)

grid.arrange(e,f,g,k,i, nrow=2, ncol=3)
grid.text("a)", x=unit(0.08, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("b)", x=unit(0.40, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("c)", x=unit(0.75, "npc"), y=unit(0.99, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("d)", x=unit(0.08, "npc"), y=unit(0.50, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
grid.text("e)", x=unit(0.40, "npc"), y=unit(0.50, "npc"), gp=gpar(fontsize=30, fontfamily="Times New Roman"))
dev.off()


---------------
  #single plots
  ---------------
  ################indv graphs
  
  png(file ="Population_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
a
dev.off()

png(file ="Seafood_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
b
dev.off()

png(file ="Landings_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
c
dev.off()

png(file ="Exploitation_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
d
dev.off()

png(file ="SST_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
e
dev.off()

png(file ="AMOPDO_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
f
dev.off()

png(file ="MEI_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
g
dev.off()

png(file ="Landings_1_v2.0_01b.png", width = 300, height = 150, units = "mm", res = 600)
h
dev.off()
png(file ="Exploitation_1_v2.0_01b.png", width = 300, height = 150, units = "mm", res = 600)
i
dev.off()
png(file ="GDP_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
j
dev.off()
png(file ="Chlorophyll_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
l
dev.off()
png(file ="Freshwater_v2.0_01.png", width = 300, height = 150, units = "mm", res = 600)
k
dev.off()
