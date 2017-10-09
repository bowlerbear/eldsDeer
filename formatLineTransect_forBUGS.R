#formatDeer_forBUGS
#formating the data from distance sampling along the line transects

#get libraries we will need
library(plyr)
library(reshape2)
library(ggplot2)

#retrieving the raw transect line data

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
datafile<-read.delim("elds_deer_LineTransectData_DBupdated.txt",header=F,as.is=T)
names(datafile)<-c("Year","Transect","TransectDist","Distance","GroupSize")
datafile$Obs<-ifelse(datafile$GroupSize==0,0,1)

#simplify transect data- ignore A and B for the moment
datafile$Transect<-gsub("A","",datafile$Transect)
datafile$Transect<-gsub("B","",datafile$Transect)
datafile$Transect<-as.numeric(datafile$Transect)

#simplify region and year
datafile$Region[datafile$Transect%in%c(1:6)]<-"Blue"
datafile$Region[datafile$Transect%in%c(7:16)]<-"Yellow"
datafile$Region[datafile$Transect%in%c(17:24)]<-"Yellow"
datafile$T<-datafile$Year-min(datafile$Year)+1

#ordering the data set
datafile<-datafile[order(datafile$Year),]
datafile<-datafile[order(datafile$Transect),]
datafile<-datafile[order(datafile$Obs,decreasing=T),]

#get summary site info
siteInfo<-unique(datafile[,c("Transect","Region")])