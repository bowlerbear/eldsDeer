#formatDeer_forBUGS
library(gdata)
library(plyr)

#read in original data file
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Div/CSVfiles")

#combile separate files in one file
allfiles<-llply(list.files(),function(x){
datafile<-read.csv(x,as.is=T,sep=";",skip=1,header=F)
datafile<-datafile[,1:11]
names(datafile)<-c("Date","Transect","TransectDistance","SightingDistance","Angle","PredPerp",
                   "Male","Female","Juvenile","Unknown","Total")

#if the whole row is empty drop it
emptyRows<-apply(datafile,1,function(x)all(x==""|is.na(x)))
datafile<-datafile[!emptyRows,]

#for Line and Length fill in blanks in a cell is empty..
for(i in 2:nrow(datafile)){
  if(datafile$Transect[i]==""){
    datafile$Transect[i]<-datafile$Transect[i-1]
  }
  if(is.na(datafile$TransectDistance[i])){
    datafile$TransectDistance[i]<-datafile$TransectDistance[i-1]
  }
}

#add year
datafile$File<-x

return(datafile)

})
datafile<-do.call(rbind,allfiles)

#Extract year
datafile$Year<-gsub("Surveys 1990-2015 line data_","",datafile$File)
datafile$Year<-as.numeric(gsub(".csv","",datafile$Year))

#simplify transect data
#datafile$Transect<-gsub("A","",datafile$Transect)#keep these in
#datafile$Transect<-gsub("B","",datafile$Transect)#keep these in
datafile$Transect<-gsub("T-","",datafile$Transect)
datafile$Transect<-gsub("\\()","",datafile$Transect)
datafile$Transect<-gsub("T","",datafile$Transect)
datafile$Transect<-gsub("\\.","",datafile$Transect)
datafile$Transect<-gsub("\\(b)","(B)",datafile$Transect)
datafile$Transect<-gsub("\\(A)","A",datafile$Transect)
datafile$Transect<-gsub("\\(B)","B",datafile$Transect)
datafile$Transect<-gsub("10 B","10B",datafile$Transect)

#recalculate angles and distances
datafile$Angle<-sapply(datafile$Angle,function(x)
  ifelse(x>180,360-x,x))
datafile$Angle[datafile$Angle==180]<-0

datafile$PerpDistance<-sin(datafile$Angle*(pi/180))*datafile$SightingDistance

#adding zeros
datafile$Male[is.na(datafile$Male)]<-0
datafile$Female[is.na(datafile$Female)]<-0
datafile$Juvenile[is.na(datafile$Juvenile)]<-0
datafile$Unknown[is.na(datafile$Unknown)]<-0

#match order of previously formatted file
datafile<-datafile[,c("Year","Transect","TransectDistance","PerpDistance","Total")]
write.table(datafile,file="elds_deer_LineTransectData_DBupdated.txt",sep="\t")

