############
#indicator calculations
#############
require(vegan)
#Pacific Islands
PI_diveristy

count<-PI_diversity[,1]
count1<-table(count)
count2<-t(count1)
count3<-as.data.frame(count2)
count3<-count3[,-1]

PIdiv<-aggregate(PI_diversity, by=PI_diversity[c("YEAR")], FUN=mean, ra.nm=T)
PIdiv<-PIdiv[,-c(1,2)]
piH<-diversity(PIdiv)
pih<-transform(piH, effective=exp(piH))
write.csv(pih, file="PI_div2.csv")

#Biomass as sum

biomass<-PI_BML1[, -c(2,4,5)]
BM<-aggregate(biomass, by=biomass[c("YEAR")], FUN=sum, ra.nm=T)
other<-PI_BML1[,-3]
O<-aggregate(other, by=other[c("YEAR")], FUN=mean, ra.nm=T)

BM<-BM[,-2]
O<-O[,-1]

PIBML<-cbind(BM, O, by="YEAR")

#All means (for mean biomass per tow)
PIBMLmeans<-aggregate(PI_BML1, by = PI_BML1[c("YEAR")], FUN=mean, ra.nm=T)

#PI drivers
PIchl<-KonIEA_CHL[! KonIEA_CHL$Chlorophyll %in% c("NaN"),]
PISeaLev<-aggregate(Kona.sea.level, by = Kona.sea.level[c("Year")], FUN=mean, ra.nm=T)
PISST<-aggregate(KonIEA_SST, by = KonIEA_SST[c("Year")], FUN=mean, ra.nm=F)
PIchla<-aggregate(PIchl, by = PIchl[c("Year")], FUN=mean, ra.nm=T)
PIenv<-cbind(PISeaLev, PISST, PIchla)
###############################
#GoMex
##############

GoMex<-merge(BTS.data, cruisedata, by="CRUISEID", all=T)
gomexL<-merge(Length.data, cruisedata, by="CRUISEID", all=T)

count2<-as.data.frame(table(gomexbio$YEAR))

gomexy<-GoMex[,-c(1,3,4,10,11,12,13,14)]

gomex<-aggregate(GoMex, by=GoMex[c("YEAR", "SPEC_BTS")], FUN = sum, na.rm=TRUE)
gomexy<-aggregate(GoMexY, by=GoMexY[c("YEAR")], FUN=mean, na.rm=TRUE)

GoMexY<-GoMexY[,-c(3,8)]
gombiom<-merge(GoMexY, count, by="CRUISE_NO.x")
gombio<-transform(gombiom, biomass=weight/Freq)
gomexbio<-aggregate(gombio, by=gombio[c("YEAR")], FUN=sum, na.rm=T)


gomexL<-merge(Length.data, cruisedata, by="CRUISEID", all=T)
GML<-subset(gomexL, MEASCD_GLF==1 & 51, select= CRUISEID:TITLE)
GML<-subset(gomexL, MEASCD_GLF==18, select= CRUISEID:TITLE)
gml<-aggregate(GML, by=GML[c("CRUISEID", "YEAR")], FUN=mean, na.rm=TRUE)
gmly<-aggregate(gml, by=gml[c("YEAR")], FUN=mean, na.rm=TRUE)

GoMexY<-aggregate(GoMex, by=GoMex[c("YEAR", "CRUISEID")], FUN = mean, na.rm=TRUE)
gomexy<-aggregate(GoMexY, by=GoMexY[c("YEAR")], FUN=mean, na.rm=TRUE)


library(dplyr)
cdtat%>%distinct(CRUISE_NO)

#################
#Biomass
####################
#Alaska
EBSBTS<-rbind(ebs201_20143, ebs1982_1984, ebs1985_1989, ebs1990_1994, ebs1995_1999, ebs2000_2004, ebs2005_2008, ebs2009_2012)
EBS<-EBSBTS[,c(4,5,7, 8,17)]
EBSt<-EBS[! EBS$NUMCPUE %in% c(-9999),]
ebs.by.haul<-aggregate(EBSt, by=EBSt[c("YEAR", "STRATUM", "HAUL")], FUN=sum, na.rm=T)
ebs.by.stratum<-aggregate(ebs.by.haul, by=ebs.by.haul[c("YEAR", "STRATUM")], FUN=mean, na.rm=T)
ebs.by.strat<-ebs.by.stratum[, c(1,2,8)]
ebs.cpue.year<-aggregate(ebs.by.strat, by = ebs.by.strat[c("YEAR")], FUN=mean, na.rm=T)
ebs.total.biomass<-transform(ebs.cpue.year, T_Biomass=WTCPUE*49901368)

#landings
ebsgf<-EBS.groundfish.2003_2013
ebsgf1<-ebsgf[, c(1,5)]
ebsgfy<-aggregate(ebsgf1, by=ebsgf1[c("YEAR")], FUN=sum, na.rm=T)

ebsgf2<-EBS.groundfish.1991_2002
ebsgf3<-ebsgf2[, c(1,5)]
ebsgfy2<-aggregate(ebsgf3, by=ebsgf3[c("YEAR")], FUN=sum, na.rm=T)

#################
#California
################

### bring in data ####
bathy00 <- data.frame(read.table("bathy.data.csv",header=TRUE,sep=",")) 
# make depth positive ####
bathy00$depth <- abs(bathy00$depth)

# set depth/lat ranges for bathymetry to get just US coast...more or less
#limit to flattery in the  north
bathy <- subset(bathy00, bathy00$lat<=48.45110 & bathy00$depth<=1200)  

# set depth lat bins
bathy$lat.bin <-      ifelse(bathy$lat < 42.5,"Blanco","Flattery")
bathy$lat.bin <-   ifelse(bathy$lat < 40.4,"Mendocino",bathy$lat.bin)
bathy$lat.bin <- 	ifelse(bathy$lat < 34.5,"Conception",bathy$lat.bin)

bathy$depth.bin <- 	ifelse(bathy$depth < 600,"slope1","slope2")
bathy$depth.bin <- 	ifelse(bathy$depth < 200,"shelf",bathy$depth.bin)

levels(as.factor(bathy$lat.bin))
levels(as.factor(bathy$depth.bin))			

# get sum of depth X lat bins
areas = aggregate(tot.area.m2 ~ depth.bin*lat.bin, FUN=sum, data=bathy)
areas$area.km2 = areas$tot.area.m2/1000000		
areas$bin = paste(areas$lat.bin,areas$depth.bin,sep="_")	

#using data from bathy5

CCBM<-aggregate(bathy5, by=bathy5[c("year")], FUN=sum, na.rm=T)

############
#NEUS
##############


#GBF<-RAF[RAF$EPU %in% c("GB"),]
#GBS<-RAS[RAS$EPU %in% c("GB"),]
#GOMF<-RAF[RAF$EPU %in% c("GOM"),]
#GOMS<-RAS[RAS$EPU %in% c("GOM"),]
#SSF<-RAF[RAF$EPU %in% c("SS"),]
#SSS<-RAS[RAS$EPU %in% c("SS"),]
#MABF<-RAF[RAF$EPU %in% c("MAB"),]
#MABS<-RAS[RAS$EPU %in% c("MAB"),]


#GBF.wide <- dcast(GBF, YEAR ~ SVSPP, value.var = "n.tow")
#GBS.wide <- dcast(GBS, YEAR ~ SVSPP, value.var = "n.tow")
#GOMF.wide <- dcast(GOMF, YEAR ~ SVSPP, value.var = "n.tow")
#GOMS.wide <- dcast(GOMS, YEAR ~ SVSPP, value.var = "n.tow")
#SSF.wide <- dcast(SSF, YEAR ~ SVSPP, value.var = "n.tow")
#SSS.wide <- dcast(SSS, YEAR ~ SVSPP, value.var = "n.tow")
#MABF.wide <- dcast(MABF, YEAR ~ SVSPP, value.var = "n.tow")
#MABS.wide <- dcast(MABS, YEAR ~ SVSPP, value.var = "n.tow")

#GBF.wide[is.na(GBF.wide)]<-0
#GBS.wide [is.na(GBS.wide)]<-0
#GOMF.wide [is.na(GOMF.wide)]<-0
#GOMS.wide [is.na(GOMS.wide)]<-0
#SSF.wide [is.na(SSF.wide)]<-0
#SSS.wide [is.na(SSS.wide)]<-0
#MABF.wide [is.na(MABF.wide)]<-0
#MABS.wide [is.na(MABS.wide)]<-0


#GBFy<-GBF.wide[,-c(1,2)]
#GBSy<-GBS.wide[,-c(1,2)]
#GOMFy<-GOMF.wide[,-c(1,2)]
#GOMSy<-GOMS.wide[,-c(1,2)]
#SSFy<-SSF.wide[,-c(1,2)]
#SSSy<-SSS.wide[,-c(1,2)]
#MABFy<-MABF.wide[,-c(1,2)]
#MABSy<-MABS.wide[,-c(1,2)]

#GBFdiv<-diversity(GBFy)
#GBSdiv<-diversity(GBSy)
#GOMFdiv<-diversity(GOMFy)
#GOMSdiv<-diversity(GOMSy)
#SSFdiv<-diversity(SSFy)
#SSSdiv<-diversity(SSSy)
#MABFdiv<-diversity(MABFy)
#MABSdiv<-diversity(MABSy)

NEFdiv<-cbind(GBFdiv, GOMFdiv, SSFdiv, MABFdiv)
NESdiv<-cbind(GBSdiv, GOMSdiv, SSSdiv, MABSdiv)


write.csv(NEFdiv, file="NEFdiv.csv")
write.csv(NESdiv, file="NESdiv.csv")

#NEsf<-cbind(NEUSFdiv, NEUSSdiv)

#NEUSdiv<-aggregate(NEsf, by=NEsf[c("YEAR")], FUN=mean, ra.nm=T)

NEUSF<-Tam_relative_abundance_fall
NEUSS<-Tam_relative_abundance_spring
NEUSraw<-rbind(NEUSF, NEUSS)

NEUSr<-NEUSraw[,-3]
NEUSR<-aggregate(NEUSr, by=NEUSr[c("YEAR", "SVSPP")], FUN=mean, ra.nm=T)
NEUSR1<-NEUSR[,-c(2,3)]
NEUS.wide <- dcast(NEUSR, YEAR ~ SVSPP, value.var = "n.tow")
NEUS.wide[is.na(NEUS.wide)]<-0
NEUSy<-NEUS.wide[,-c(1)]
NEUSdiv<-diversity(NEUSy)
NEUSeff<-transform(NEUSdiv, effective=exp(NEUSdiv))
write.csv(NEUSeff, file="NEUSdiversity.csv")

RAF1<-RAF[,-3]
RAS1<-RAS[,-3]

NEUSF<-aggregate(RAF1, by=RAF1[c("YEAR", "SVSPP")], FUN=sum, ra.nm=F)
NEUSS<-aggregate(RAS1, by=RAS1[c("YEAR", "SVSPP")], FUN=sum, ra.nm=F)

NEUSF.wide <- dcast(NEUSF, YEAR ~ SVSPP, value.var = "n.tow")
NEUSS.wide <- dcast(NEUSS, YEAR ~ SVSPP, value.var = "n.tow")

NEUSF.wide[is.na(NEUSF.wide)]<-0
NEUSS.wide[is.na(NEUSS.wide)]<-0

NEUSFy<-NEUSF.wide[,-c(1,2)]
NEUSSy<-NEUSS.wide[,-c(1,2)]



GBFef<-transform(GBFdiv, effective=exp(GBFdiv))
GBSef<-transform(GBSdiv, effective=exp(GBSdiv))
GOMFef<-transform(GOMFdiv, effective=exp(GOMFdiv))
GOMSef<-transform(GOMSdiv, effective=exp(GOMSdiv))
SSFef<-transform(SSFdiv, effective=exp(SSFdiv))
SSSef<-transform(SSSdiv, effective=exp(SSSdiv))
MABFef<-transform(MABFdiv, effective=exp(MABFdiv))
MABSef<-transform(MABSdiv, effective=exp(MABSdiv))

NEFef<-cbind(GBFef, GOMFef, SSFef, MABFef)
NESef<-cbind(GBSef, GOMSef, SSSef, MABSef)

write.csv(NEFef, file="NEFef.csv")
write.csv(NESef, file="NESef.csv")

count1<-table(NEUSR)

#NEUS<-diversity(PIdiv)
#pih<-transform(piH, effective=exp(piH))
#write.csv(pih, file="PI_div2.csv")
