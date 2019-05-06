NEUSbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/NEUSbio.csv")
NEUSenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/NEUSenv.csv")
EBSbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/EBSbio.csv")
EBSenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/EBSenv.csv")
ccbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/ccbio.csv")
ccenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/ccenv.csv")
GOMbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/GOMbio.csv")
GOMenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/GOMenv.csv")
PIbio <- read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/PIbio.csv")
PIenv<-read.csv("~/jtwork/IEAnatThresh/GFResults2016_02/PIenv.csv")

PIenvNO12<-PIenv[-12,]

NEUSdat<-cbind(NEUSbio, NEUSenv)
EBSdat<-cbind(EBSbio, EBSenv)
ccdat<-cbind(ccbio, ccenv)
GOMdat<-cbind(GOMbio, GOMenv)
PIdat<-cbind(PIbio, PIenvNO12)

GDP<-cbind(NEUSdat$GDP, EBSdat$GDP, ccdat$GDP, PIdat$GDP)

plot(NEUSdat$GDP.inc)
plot(ccdat$Rich~ccdat$Freshwater)

ccglm1<-glm(ccdat$Length~ccdat$Freshwater)