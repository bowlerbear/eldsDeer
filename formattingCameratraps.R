##################################################################################################
#
#file to retreive and begin formatting the observations of deer in the camera trapping
#
##################################################################################################

#libraries
library(maptools)
library(rgdal)
library(raster)
library(plyr)
library(ggplot2)
library(lubridate)
library(reshape2)

###################################################################################################

#get shape file of the Sanctuary
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/MMR layers/ProtectedArea")
sws<-readOGR(dsn=getwd(),layer="WDPA_Aug2017_protected_area_1236-shapefile-polygons")

#get location of transects
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/MMR layers/Updated")
transects<-readOGR(getwd(),layer="Transects")

#give them the transect names
plot(transects)
transects$Id<-as.character(transects$Id)
transects$Id[1:16]<-1:16
transects$Id[17:24]<-24:17

######################################################################################################

setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/Camera traps")

###############
#station infos#
###############

#stations 2014
stations2014<-read.delim("Data 2014 & 2016 Clean 011216_2014stations.txt",skip=1,as.is=T,dec=",")
stations2014Locations<-unique(stations2014[,c("GRID.CELL.ID","CAMERA.ID","SD.CARD.ID.SETUP","LATITUDE..N.",
                                              "LONGITUDE..E.","LURE.TYPE","WATER.POINT","SALT.LICK",
                                              "HABITAT","TERRAIN","STATION.POTENTIAL","SETUP.DATE",
                                              "TAKE.DOWN.DATE","PATH.TYPE","CAMERA.OPERATIONAL.TAKE.DOWN")])
stations2014Locations$Year<-2014

#stations2016
stations2016<-read.delim("Data 2014 & 2016 Clean 011216_2016stations.txt",skip=1,as.is=T,dec=",")
stations2016Locations<-unique(stations2016[,c("GRID.CELL.ID","CAMERA.ID","SD.CARD.ID.SETUP","LATITUDE..N.",
                                              "LONGITUDE..E.","LURE.TYPE","WATER.POINT","SALT.LICK",
                                              "HABITAT","TERRAIN","STATION.POTENTIAL","SETUP.DATE",
                                              "TAKE.DOWN.DATE","PATH.TYPE","CAMERA.OPERATIONAL.TAKE.DOWN")])
stations2016Locations$Year<-2016

#all stations
locations<-rbind(stations2014Locations,stations2016Locations)
locations<-subset(locations,!is.na(LATITUDE..N.)|!is.na(LONGITUDE..E.))

###########################################################
#get forest cover and elevational data at all lat and lons#
###########################################################
locationsSP<-locations
coordinates(locationsSP)<-c("LONGITUDE..E.","LATITUDE..N.")
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/MMR layers/Updated")
top<-raster("MMR_alt.grd")
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/MMR layers/Updated/Hansen")
forest2000<-raster("Hansen_GFC2015_treecover2000_30N_090E.tif")
locations$ForestCover<-extract(forest2000,locationsSP)
locations$Topography<-extract(top,locationsSP)

##########################
#TIDY THE STATION INFO UP#
##########################

#remove "SD" from SD.CARD.ID.SETUP
locations$SD.CARD.ID.SETUP<-as.numeric(gsub("SD","",locations$SD.CARD.ID.SETUP))
locations$CAMERA.ID<-gsub("_","",locations$CAMERA.ID)
locations$CAMERA.ID<-gsub("_","",locations$CAMERA.ID)

#format dates
locations$SETUP.DATE<-gsub("des","dec",locations$SETUP.DATE)
locations$SETUP.DATE<-as.Date(locations$SETUP.DATE,format="%d-%b-%y")
locations$Month<-month(locations$SETUP.DATE)
locations$TAKE.DOWN.DATE<-as.Date(locations$TAKE.DOWN.DATE,format="%d-%b-%y")
locations$TAKE.DOWN.DATE[which(locations$TAKE.DOWN.DATE=="2015-02-06")]<-"2016-02-06"#fix the year

#note: there isnt always a take down date...

#clean the variables into standard classes
locations$LURE.TYPE[locations$LURE.TYPE%in%c("CASTOR","CASTOR + FISH OIL","CATNIP")]<-"Yes"
locations$WATER.POINT[locations$WATER.POINT%in%c("FALSE","No")]<-"No"
locations$WATER.POINT[locations$WATER.POINT%in%c("TRUE","YES")]<-"Yes"
locations$SALT.LICK[locations$SALT.LICK%in%c("FALSE","No","no")]<-"No"
locations$SALT.LICK[locations$SALT.LICK%in%c("TRUE")]<-"Yes"
locations$SALT.LICK[locations$SALT.LICK%in%c("?")]<-NA
locations$STATION.POTENTIAL[locations$STATION.POTENTIAL%in%c("BAD","POOR")]<-"POOR"
locations$HABITAT[locations$HABITAT%in%c("INDAING","INDAING FOREST")]<-"INDAING FOREST"
locations$HABITAT[locations$HABITAT%in%c("Upper mixed deciduous","UPPER MIXED DECIDUOUS",
                                 "UPPER MIXED DECIDUOUS FOREST")]<-"UPPER MIXED DECIDUOUS FOREST"
locations$HABITAT[locations$HABITAT%in%c("DRY FOREST","EDRY FOREST")]<-"DRY FOREST"

#######################
#get deer observations#
#######################

#reading in each file on the camera trap observations
#and subsetting to only the elds deer observations
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/Camera traps")
deer2014<-read.delim("Data 2014 & 2016 Clean 011216_2014.txt",as.is=T)
deer2014<-subset(deer2014,Wild.Animals%in%c("Rucervus eldii"))#142
deer2016<-read.delim("Data 2014 & 2016 Clean 011216_2016.txt",as.is=T)
deer2016<-subset(deer2016,Wild.Animals%in%c("Rucervus eldii"))#35
deer2016_2<-read.delim("Data 2014 & 2016 Clean 011216_2016_2.txt",as.is=T)
deer2016_2<-subset(deer2016_2,Wild.Animals%in%c("Rucervus eldii"))#112

#combining each file for deer
myvars<-c("Year","Location","Camera","SD.Card","X.10","Start...Event","End...Event")
deer2014<-deer2014[,myvars]
deer2016<-deer2016[,myvars]
deer2016_2<-deer2016_2[,myvars]
deer<-rbind(deer2014,deer2016,deer2016_2)

#tidying the file names a little
deer$X.10<-as.numeric(deer$X.10)
names(deer)[5]<-"numberDeer"
names(deer)[6:7]<-c("StartDate","EndDate")
deer$StartDate<-as.Date(deer$StartDate,format="%d.%m.%Y")
deer$Year[which(deer$Year==2015)]<-2016#put the december 2015 points into 2016 season

#formatting the variables so they the columns can be merged with the station location data 
deer$Location<-gsub("_","",deer$Location)
deer$Camera<-gsub("_","",deer$Camera)
deer$Location<-gsub("-","",deer$Location)
deer$Location<-gsub("A","",deer$Location)
deer$Camera<-gsub("-","",deer$Camera)
deer$Location<-gsub("A","",deer$Location)
deer$Location<-gsub("b","",deer$Location)
locations$GRID.CELL.ID<-gsub("A","",locations$GRID.CELL.ID)
locations$GRID.CELL.ID<-gsub("B","",locations$GRID.CELL.ID)
locations$GRID.CELL.ID[which(locations$GRID.CELL.ID=="In military area, outside grid")]<-"SWSRMY"
locations$GRID.CELL.ID<-gsub(" ","",locations$GRID.CELL.ID)
locations$CAMERA.ID<-gsub("NINA0","NINA",locations$CAMERA.ID)
locations$CAMERA.ID<-gsub("nina","NINA",locations$CAMERA.ID)

#create ID for merging the deer observation data with the camera trap station data
deer$ID<-paste(deer$Year,deer$Location,deer$Camera,deer$SD.Card,sep="_")
locations$ID<-paste(locations$Year,locations$GRID.CELL.ID,locations$CAMERA.ID,locations$SD.CARD.ID.SETUP,sep="_")

#still a few grids for which we have observations are not recorded as a station??
unique(deer$ID[!deer$ID%in%locations$ID])
#[1] "2016_SWS104_NINA22_26" "2016_SWS105_NINA15_11"
#[3] "2016_SWS76_NINA13_14"
#We loose 8 observations of 1 or 2 deer

#table(deer$Location)
#majority of deer records came from grid cell 119 

##############################################
#create a list of survey days for each camera#
##############################################

surveys<-locations[,c("ID","SETUP.DATE","TAKE.DOWN.DATE","Year")]
#surveys$SurveyDays<-surveys$TAKE.DOWN.DATE-surveys$SETUP.DATE
#tapply(surveys$SurveyDays,surveys$Year,median,na.rm=T)
#2014 2016 
#10   38

#assume missing surveys in 2016 were for 38 days
surveys$TAKE.DOWN.DATE[is.na(surveys$TAKE.DOWN.DATE)&surveys$Year==2016]<-surveys$SETUP.DATE[is.na(surveys$TAKE.DOWN.DATE)&surveys$Year==2016]+38
surveys$TAKE.DOWN.DATE[is.na(surveys$TAKE.DOWN.DATE)&surveys$Year==2014]<-surveys$SETUP.DATE[is.na(surveys$TAKE.DOWN.DATE)&surveys$Year==2014]+10
surveys$nDays<-as.numeric(surveys$TAKE.DOWN.DATE-surveys$SETUP.DATE)

#for each survey date add the day number
library(plyr)
surveys<-ddply(surveys,.(ID),function(x){
  Day=1:max(x$nDays)
  data.frame(Day=Day,x)
})
surveys$Date<-surveys$SETUP.DATE+surveys$Day-1
#the surveys dates can now be merged with the dates when deer were seen

#aggregate the deer information to numbers seen per day
deer<-ddply(deer,.(ID,StartDate),summarise,numberDeer=sum(numberDeer))

#create ids in both the deer and survey data sets so they can be matching
deer$StartDateID<-paste(deer$StartDate,deer$ID,sep="_")
surveys$StartDateID<-paste(surveys$Date,surveys$ID,sep="_")

#Add the deer observations to the survey/location info dataset
surveys$numberDeer<-deer$numberDee[match(surveys$StartDateID,deer$StartDateID)]
surveys$numberDeer[is.na(surveys$numberDeer)]<-0
surveys$PA<-ifelse(surveys$numberDeer>0,1,0)

table(surveys$PA)
#0    1 
#4258  106

table(surveys$numberDeer)
#mostly 0 or 1

#get list of ID with potential operational issues
probIDs<-locations$ID[locations$CAMERA.OPERATIONAL.TAKE.DOWN%in%c("No","FALSE","No (Lost")]
#keep in one from that cell
#unique(locations[,c("ID","LONGITUDE..E.","LATITUDE..N.")])
probIDs<-probIDs[-6]
surveys2<-subset(surveys,ID%in%probIDs)
sum(surveys2$PA)#nothing was seen at these stations, so remove them
surveys<-subset(surveys,!ID%in%probIDs)
#####################################################################################

