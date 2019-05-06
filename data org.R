##organizing data from OceanAdapt (NEUS)

require(reshape2)
require(XLConnect)
require(ggplot2)
require(plyr)
require(vegan)
require(RColorBrewer)
require(grid)

#Removes unwanted columns
neusdat<-neus_data[,-c(1,2,3,5,6,9,10,11,12,13,14,15)]

#sort and aggregate data
ALland<-aggregate (Alaska.landings, by = Alaska.landings[c("year", "area")], FUN =mean, na.rm=TRUE)
RDyear<-aggregate(richdivcount, by = richdivcount[c("year")], FUN = mean, na.rm=TRUE)
richyear2003<-richdivcount[richdivcount$year=="2003",]

#Diversity
div<-diversity(RDyear, index = "simpson")
rich<-rich(matrix = richyear, nrandom=499, verbose = TRUE)