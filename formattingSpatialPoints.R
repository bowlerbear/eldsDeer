#compile a list of the folders containing the data each year

#set working directory containing a folder for each year
base<-"C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Reports/Copies of SWS Census reports/Copies of SWS Census reports"
setwd(base)
myfolders<-list.files()


#check that each year has a georeferenced folder
lapply(myfolders,function(x)any(grepl("georeferenced",list.files(paste(getwd(),x,sep="/")))))
# 2003 onwards we do

#restrict folders to 2003 onwards
myfolders<-myfolders[8:length(myfolders)]

#now check that each of these folders has a spatial points shape file
lapply(myfolders,function(x)any(grepl("SpatialPoints.shp",list.files(paste(getwd(),x,"georeferenced",sep="/")))))
#yep!!

#open up each folder and get the spatial points objects
library(rgdal)
library(plyr)
out<-ldply(myfolders,function(x){
  mydirectory<-paste(base,x,"georeferenced",sep="/")
  mypoints<-readOGR(dsn=mydirectory,layer="SpatialPoints")
  mydata<-data.frame(mypoints@coords)
  mydata$Year<-as.numeric(x)
  return(mydata)
})


#writing the file
write.table(out,file="SpatialPoints_BasicOutput.txt",sep="\t")

#add the counts off the maps

#read in procesed file
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Reports")
mydata<-read.delim("SpatialPoints_ProcessedOutput.txt")
library(ggplot2)
qplot(coords.x1,coords.x2,data=mydata)+
  facet_wrap(~Year)+
  geom_point(aes(colour=log(Count)),size=5)+
  scale_colour_continuous(low="lightgray",high="red")+
  theme_bw()


#looking at the raw data file
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
datafile<-read.delim("elds_deer_LineTransectData.txt",header=F,as.is=T)
names(datafile)<-c("Year","Transect","TransectDist","Distance","GroupSize","Region")
datafile$GroupSize[is.na(datafile$GroupSize)]<-0
datafile$Obs<-ifelse(datafile$GroupSize==0,0,1)


#plotting total seen per year
library(plyr)
totalCounts<-ddply(datafile,.(Year,Region),summarize,tot=sum(GroupSize))
library(ggplot2)
qplot(Year,tot,data=totalCounts,colour=Region)

#################################################################################