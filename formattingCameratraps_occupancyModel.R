#get the file that already has been partly formatted

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/GitFiles/formattingCameratraps.R')

#########################################################
# using Grid cell as the site# 
##########################################################

#pull out the grid and year data from the ID variable
surveys$Grid<-sapply(surveys$ID,function(x)strsplit(x,"_")[[1]][2])
surveys$Year<-sapply(surveys$ID,function(x)strsplit(x,"_")[[1]][1])
surveys$Camera<-sapply(surveys$ID,function(x)strsplit(x,"_")[[1]][3])

#recount Day so it goes over multiple years/cameras in the same grid cells
surveys<-ddply(surveys,.(Grid),function(x){
  Day=1:nrow(x)
  data.frame(cbind(Day2=Day,x))
})

#get rid of data without water information
nowaterData<-locations$ID[locations$WATER.POINT==""]
surveys<-subset(surveys,!ID%in%nowaterData)

#convert to an array
obsMatrix<-acast(surveys,Grid~Day2,value.var="PA")

#Restrict to the first set of rows
nu<-30
y<-obsMatrix[,1:nu]
#nrow(y[apply(y,1,function(x)any(x)!=0),])
#51 site/year observations..

#get covariate data

#occupancy factors
occupancyCovariates<-unique(locations[,c("GRID.CELL.ID","ForestCover")])
occupancyCovariates<-ddply(occupancyCovariates,.(GRID.CELL.ID),summarise,ForestCover=mean(ForestCover))
occupancyCovariates$ForestCoverSQD<-(occupancyCovariates$ForestCover+0.001)^2
#occupancyCovariates[,2:3]<-sapply(occupancyCovariates[,2:3],scale)

#detection covariates
detectionCovariates<-unique(locations[,c("ID","GRID.CELL.ID","Year","WATER.POINT")])
detectionCovariates<-merge(surveys,detectionCovariates,by=c("ID","Year"),all.x=T)
detectionCovariates$Year<-as.numeric(detectionCovariates$Year)
detectionCovariates$WATER.POINT<-ifelse(detectionCovariates$WATER.POINT=="Yes",1,0)
WATER<-acast(detectionCovariates,Grid~Day2,value.var="WATER.POINT")[,1:nu]
YEAR<-acast(detectionCovariates,Grid~Day2,value.var="Year")[,1:nu]

#put any value for the missing data
WATER[is.na(WATER)]<-0
YEAR[is.na(YEAR)]<-2016
YEAR[YEAR==2014]<-0
YEAR[YEAR==2016]<-1

#make sure we have all data
nrow(occupancyCovariates)==nrow(y)
occupancyCovariates<-occupancyCovariates[order(match(occupancyCovariates$GRID.CELL.ID,row.names(y))),]
WATER<-WATER[order(match(row.names(WATER),row.names(y))),]
YEAR<-YEAR[order(match(row.names(YEAR),row.names(y))),]
