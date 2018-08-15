#################################################################################

#run file to retreive and begin formatting the observations of deer in the camera trapping
source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingCameratraps.R')

#################################################################################

#get additional libraries we will need
library(mgcv)
library(RColorBrewer)

#################################################################################

#Making the grid

#import 2 X 2 km grid giving the SWS grid ids and make it into a raster

#the original grid of the study
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Camera traps")
grid<-read.delim("CameraGrid.txt")
grid<-subset(grid,Col!=19)

#make new grid
gridExtent<-c(94.35,94.69,19.98,20.36)
rO<-raster(extent(gridExtent),nrow=21,ncol=18)#divisable by 3
rO[]<-1:ncell(rO)
projection(rO)<-CRS("+proj=longlat +datum=WGS84")
#convert it into a dataframe
rODF<-as.data.frame(rO,xy=T)
names(rODF)[which(names(rODF)=="layer")]<-"cellnu"
#also add on grid info
rODF$swsGrid<-grid$Grid

#map this to our new grid

#the larger 3km grid - approximating the home range size of the deer
temp<-projectExtent(rO,crs=CRS("+proj=eck4"))
extent(temp)
#xmin        : 8646427 
#xmax        : 8685474 
#ymin        : 2601521 
#ymax        : 2649731
myGrid<-raster(extent(c(8655000,8680000,2613000,2645000)),crs=CRS("+proj=eck4"))
res(myGrid)<-3000
myGrid[]<-1
myGrid<-projectRaster(myGrid,crs="+proj=longlat")
myGrid[]<-1:ncell(myGrid)
myGrid3km<-myGrid

#plotting
plot(myGrid)
plot(sws,add=T)

#and then to the smaller 1km grid  - well within the distance caught by the camera traps and transects
myGrid<-disaggregate(myGrid,fact=3)
myGridDF<-as.data.frame(myGrid,xy=T)
#now get its smaller numbers
myGridDF$Grid1km<-1:ncell(myGrid)
names(myGridDF)[3]<-"Grid3km"
myGrid[]<-1:ncell(myGrid)#have the smaller 1km grid as the layer cells

#plotting
plot(myGrid)
plot(sws,add=T)
plot(transects,add=T)

#for later use, lets call myGrid as r
r<-myGrid

#get list of grid cells we want to impute for - covered by at least 50% of the reserve
myGridMask<-data.frame(extract(myGrid3km,sws,weights=T,normalizeWeights=F)[[1]])
myGridMask<-myGridMask$value[myGridMask$weight>0.55]
myGrid2<-myGrid3km
myGrid2[!getValues(myGrid3km)%in%myGridMask]<-NA
plot(myGrid2)
plot(sws,add=T)
#looks good!
grid3kmImpute <- getValues(myGrid2)[!is.na(getValues(myGrid2))]
  
#############################################################
#Get species presence data from camera traps or line methods#
#############################################################

#get spatial points
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Reports")

#read in spatial point data with distances
transectPoints<-read.delim("SpatialPoints_ProcessedOutput(2013-2017distances).txt")
transectPoints<-subset(transectPoints,Year%in%c(2014:2016))
coordinates(transectPoints)<-c("coords.x1","coords.x2")

#get camera trap points
#standard occupancy analysis on presence/absence at area surrouding the camera traps
#run the first section of code in 'formating_Cameratraps_for_BUGS.R'
surveysP<-subset(surveys,numberDeer>0)
surveysP<-merge(surveysP,locations,by="ID")
surveysP$Site<-paste(surveysP$LONGITUDE..E.,surveysP$LATITUDE..N.)
coordinates(surveysP)<-c("LONGITUDE..E.","LATITUDE..N.") 

#Identify cells with a presence observation
r[]<-1:ncell(r)
transectPointsRP<-extract(r,transectPoints)
surveysRP<-extract(r,surveysP)
r[]<-0
r[transectPointsRP]<-1
r[surveysRP]<-1
plot(r)
plot(sws,add=T)

######################
#Get absence data too#
#####################

#from the camera traps
locations<-subset(locations,!ID%in%probIDs)#remove those without functioning camera on collection (and never recorded anything)
locations$Site<-paste(locations$LONGITUDE..E.,locations$LATITUDE..N.)
surveysRA<-subset(locations,!Site%in%surveysP@data$Site)
coordinates(surveysRA)<-c("LONGITUDE..E.","LATITUDE..N.")
proj4string(surveysRA)<-crs(crs(sws))

#remove those beyond the edge of the boundary
library(rgeos)
surveysRA<-gIntersection(surveysRA,sws)
plot(sws)
plot(surveysRA,add=T)

#identify the cells
r[]<-1:ncell(r)
surveysRA<-extract(r,surveysRA)

#from the transects
transectsRA<-do.call(c,extract(r,transects))
transectsRA<-transectsRA[!transectsRA%in%transectPointsRP]

par(mfrow=c(1,3))
#camera trap alone
r[]<-0
r[surveysRA]<-1
r[surveysRP]<-4
plot(r,col=c("white","azure3","red"),axes=FALSE,legend=FALSE,box=FALSE)
plot(sws,add=T)

#transect alone
r[]<-0
r[transectsRA]<-2
r[transectPointsRP]<-5
plot(r,col=c("white","azure3","red"),axes=FALSE,legend=FALSE,box=FALSE)
plot(sws,add=T)

#plotting presences and absences together
r[]<-0
r[surveysRA]<-1
r[transectsRA]<-2
r[surveysRP]<-4
r[transectPointsRP]<-5
plot(r,col=c("white","azure3","azure3","red","red"),axes=FALSE,legend=FALSE,box=FALSE)
plot(sws,add=T)

#######################
#Obtain the covariates#
#######################

#group together the grid cell observations into transect lines

#is a cell inside or outside the military area
r[]<-1:ncell(r)
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers")
military<-readOGR(getwd(),layer="Military area")
militaryR<-ldply(extract(r,military,weights=T,normalizeWeights=F))

#get average forest area for each cell
fdir<-"C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated/Hansen"
forest2000<-raster(paste(fdir,"Hansen_GFC2015_treecover2000_30N_090E.tif",sep="/"))
forest2000<-crop(forest2000,extent(r))
forest2000M<-mask(forest2000,sws)
#plot(forest2000M,col=brewer.pal(9,"Greys"))
#plot(sws,add=T)
forest2000<-as.data.frame(forest2000,xy=T)
coordinates(forest2000)<-c("x","y")
forest2000@data$Cell<-extract(r,forest2000)
#get average cover per cell
forestDF<-forest2000@data
forestDF<-ddply(forestDF,.(Cell),summarise,ForestCover=mean(Hansen_GFC2015_treecover2000_30N_090E,na.rm=T))
#Tree cover in the year 2000, defined as canopy closure for all vegetation taller than 5m in height. 
#Encoded as a percentage per output grid cell, in the range 0–100.

#import Johns work in google earth
fields<-readOGR(getwd(),layer="JohnsPolygons")
fields<-subset(fields,grepl("Field",fields$Name))
fields<-subset(fields,!grepl("military",fields$Name))
library(rgeos)
fields<-gUnaryUnion(fields)
#plot(sws)
#plot(fields,add=T,col="grey")
#calculate coverage of each grid cell by a field
r[]<-1:ncell(r)
fieldsR<-ldply(extract(r,fields,weights=T,normalizeWeights=F))
fieldsR<-ddply(fieldsR,.(value),summarise,weight=sum(weight))
#fraction of cell covering by fields

#also get information on nearby villages
villages<-readOGR(getwd(),layer="JohnsPolygons")
villages<-subset(villages,grepl("Village",villages$Name))
villages<-gUnaryUnion(villages)
plot(sws)
plot(villages,add=T,col="grey")
#extract density
#probably should increase the extent of r
r2<-r
extent(r2)<-c((extent(r)[1]-(10*res(r)[1])),(extent(r)[2]+(10*res(r)[1])),
              (extent(r)[3]-(10*res(r)[2])),(extent(r)[4]+(10*res(r)[2])))
res(r2)<-res(r)
r2[]<-1:ncell(r2)
villagesR<-ldply(extract(r2,villages,weights=T,normalizeWeights=F))
#sum it up over each cell
villagesR<-ddply(villagesR,.(value),summarise,weight=sum(weight))

#overlay onto the raster
rVillages<-as.data.frame(r2,xy=T)
rVillages$Village<-villagesR$weight[match(rVillages$layer,villagesR$value)]
rVillages$Village[is.na(rVillages$Village)]<-0
villageRaster<-rasterFromXYZ(rVillages[,c("x","y","Village")])##fraction to which a cell is covered by a village

#smooth the raster
villageRaster <- focal(villageRaster, w=matrix(1, 15, 15), mean,na.rm=T)
projection(villageRaster)<-CRS("+proj=longlat +ellps=WGS84") 
villageRaster<-projectRaster(from=villageRaster,to=myGrid)
rVillages<-as.data.frame(r2,xy=T)
plot(villageRaster)
plot(sws,add=T)
plot(villages,add=T)

#convert to data frame
villagesDF<-as.data.frame(villageRaster)
villagesDF$cell<-1:ncell(villageRaster)

##################
#combine all data#
##################

#on the 1km grid

#myGridDF$Military<-sapply(myGridDF$Grid1km,function(x)ifelse(x%in%militaryR,1,0))
myGridDF$Military<-militaryR$weight[match(myGridDF$Grid1km,militaryR$value)]
myGridDF$ForestCover<-forestDF$ForestCover[match(myGridDF$Grid1km,forestDF$Cell)]
myGridDF$ForestCover<-myGridDF$ForestCover/100
myGridDF$Fields<-fieldsR$weight[match(myGridDF$Grid1km,fieldsR$value)]
myGridDF$Fields[is.na(myGridDF$Fields)]<-0
myGridDF$Villages<-villagesDF$layer[match(myGridDF$Grid1km,villagesDF$cell)]
#check<-rasterFromXYZ(myGridDF[,c("x","y","Villages")])

#add deer data from camera traps
r[]<-NA
r[c(surveysRA)]<-0
r[c(surveysRP)]<-1
DF<-as.data.frame(r)
myGridDF$Deer.ct<-DF$layer

#add deer data from transects
r[]<-NA
r[c(transectsRA)]<-0
r[c(transectPointsRP)]<-1
DF<-as.data.frame(r)
myGridDF$Deer.lt<-DF$layer

#identify squares with a transect
r[]<-NA
r[transectsRA]<-1
DF<-as.data.frame(r)
myGridDF$Transect<-DF$layer

#also organise the data on a 3km grid
myGridDF3km<-ddply(myGridDF,.(Grid3km),summarise,
                   Military=mean(Military,na.rm=T),
                   Villages=mean(Villages,na.rm=T),
                   ForestCover=mean(ForestCover,na.rm=T),
                   Fields=mean(Fields,na.rm=T),
                   Deer.ct=max(Deer.ct,na.rm=T),
                   Deer.lt=max(Deer.lt,na.rm=T),
                   TransectNu=sum(Transect,na.rm=T))
#add NAs
myGridDF3km$Military[is.na(myGridDF3km$Military)]<-0
myGridDF3km$Deer.ct[is.infinite(myGridDF3km$Deer.ct)]<-NA
myGridDF3km$Deer.lt[is.infinite(myGridDF3km$Deer.lt)]<-NA

#add a general deer presence value
myGridDF3km$Deer<-ifelse(myGridDF3km$Deer.ct==1|myGridDF3km$Deer.lt==1,1,0)
#if one method says 0, then lets say it is zero
myGridDF3km$Deer[is.na(myGridDF3km$Deer)]<-0

#making factors
myGridDF3km$ForestCoverF<-ifelse(myGridDF3km$ForestCover<0.01,0,1)
myGridDF3km$ForestCoverF2<-ifelse(myGridDF3km$ForestCover<0.40,0,1)
myGridDF3km$MilitaryF<-ifelse(myGridDF3km$Military>0,1,0)

########################################################################

#Fitting a hierarcical model to both datasets

#get full list of grid cells with data

##########################
#(1) get camera trap data#
##########################

#initally use the 1km grid
cameraTraps<-merge(surveys,locations,by="ID")
coordinates(cameraTraps)<-c("LONGITUDE..E.","LATITUDE..N.")
proj4string(cameraTraps)<-crs(crs(sws))
cameraTraps<-cameraTraps[!is.na(over(as(cameraTraps,'SpatialPoints'),as(sws,'SpatialPolygons'))),]
cameraTraps$Cell<-extract(myGrid,cameraTraps)
cameraTraps<-merge(myGridDF,data.frame(cameraTraps),by.x="Grid1km",by.y="Cell")
cameraTraps$Year<-sapply(as.character(cameraTraps$ID),function(x)strsplit(x,"_")[[1]][1])

#get rid of data without water information
#nowaterData<-locations$ID[locations$WATER.POINT==""]
#cameraTraps<-subset(cameraTraps,!ID%in%nowaterData)

#combine water point or salt lick as an attractor
cameraTraps$LURE.TYPE<-as.character(cameraTraps$LURE.TYPE)
cameraTraps$SALT.LICK<-as.character(cameraTraps$SALT.LICK)
cameraTraps$WATER.POINT<-as.character(cameraTraps$WATER.POINT)
cameraTraps$LURE.TYPE[cameraTraps$LURE.TYPE=="?"]<-"NONE"
cameraTraps$LURE.TYPE[cameraTraps$LURE.TYPE==""]<-"NONE"
cameraTraps$WATER.POINT[cameraTraps$WATER.POINT==""]<-"No"
cameraTraps$SALT.LICK[cameraTraps$SALT.LICK==""]<-"No"
attractorDF<-unique(cameraTraps[,c("ID","Grid1km","LURE.TYPE","WATER.POINT","SALT.LICK","Year")])

#how many stations have attractors
table(attractorDF$LURE.TYPE)
table(attractorDF$WATER.POINT)
table(attractorDF$SALT.LICK)#only 8

#also add a yes if there is a salt lick
cameraTraps$WATER.POINT[cameraTraps$SALT.LICK=="Yes"]<-"Yes"

#check attractor is constant within a grid cell
out<-ddply(attractorDF,.(Grid1km),summarise,nu=length(unique(WATER.POINT)))
#5 still have variation in whether it is a water point
dups<-subset(attractorDF,Grid1km%in%out$Grid1km[out$nu==2])
#take those from 2016
cameraTraps<-subset(cameraTraps,!cameraTraps$ID%in%dups$ID[dups$Year=="2014"])

#check lure is constant with a cell
out<-ddply(attractorDF,.(Grid1km),summarise,nu=length(unique(LURE.TYPE)))
dups<-subset(attractorDF,Grid1km%in%out$Grid1km[out$nu==2])
#take those from 2016
cameraTraps<-subset(cameraTraps,!cameraTraps$ID%in%dups$ID[dups$Year=="2014"])

#in each 1km grid - how many sites were there
#ddply(cameraTraps,.(Grid1km),summarise,nu=length(unique(ID))) - usually just 1

#plot presences
presDF<-ddply(cameraTraps,.(x,y),summarise,P=max(PA,na.rm=T))
presDF<-subset(presDF,P==1)
coordinates(presDF)<-c("x","y")
plot(sws)
plot(presDF,add=T)

#recount Day so it goes over multiple years/cameras in the same grid cells
cameraTraps<-ddply(cameraTraps,.(Grid1km),function(x){
  x<-x[order(x$Year,decreasing = T),]
  Day=1:nrow(x)
  data.frame(cbind(Day2=Day,x))
})

#relabel Grid to be 1 to 9 within each Grid3km
cameraTraps<-ddply(cameraTraps,.(Grid3km),function(x){
  x$GridRep=as.factor(x$Grid1km)
  levels(x$GridRep)<-1:length(unique(x$GridRep))
  x$GridRep<-as.numeric(x$GridRep)
  return(x)
})

#only use data from our cells of interest
cameraTraps<-subset(cameraTraps,Grid3km%in%grid3kmImpute)

#convert data into a 3-dimensional array
obsMatrix<-acast(cameraTraps,Grid3km~GridRep~Day2,value.var="PA")
dim(obsMatrix)

nu<-30
y<-obsMatrix[,,1:nu]

#water info
cameraTraps$WATER.POINT<-sapply(cameraTraps$WATER.POINT,function(x)ifelse(x=="Yes",1,0))
waterMatrix<-acast(cameraTraps,Grid3km~GridRep,value.var="WATER.POINT",fun=function(x) 
  round(mean(x)))
waterMatrix[is.na(waterMatrix)]<-0

#lure info
cameraTraps$LURE.TYPE<-sapply(cameraTraps$LURE.TYPE,function(x)ifelse(x=="Yes",1,0))
lureMatrix<-acast(cameraTraps,Grid3km~GridRep,value.var="LURE.TYPE",fun=function(x)
  round(mean(x)))
lureMatrix[is.na(lureMatrix)]<-0

#month info
cameraTraps$Month<-sapply(cameraTraps$Month,function(x)ifelse(x==12,0,x))
monthMatrix<-acast(cameraTraps,Grid3km~GridRep~Day2,value.var="Month",fun=max)[,,1:nu]
monthMatrix[is.infinite(monthMatrix)]<-0

#list of camera trap sites
sites.ct<-acast(cameraTraps,Grid3km~GridRep~Day2,value.var="PA")
sites.ct<-as.numeric(dimnames(sites.ct)[[1]])

###########################################
#(2) get spatial points along the transect#
###########################################

#get presences
transectPoints$cells<-extract(myGrid3km,transectPoints)
transectPoints<-transectPoints[,c("Count","Distance","cells","Year")]
plot(sws)
plot(transectPoints,add=T,col="red")

#convert into a data frame
transectDF<-transectPoints@data
transectDF$GroupSize<-transectDF$Count
#truncate observation at 250 m
transectDF<-subset(transectDF,Distance<=250)
nrow(transectDF)#76 detection events

#get transect effort

#read in transects for each year
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Reports/Copies of SWS Census reports/Copies of SWS Census reports/2016/georeferenced"
transect2016<-readOGR(tdir,layer="TransectLines2016")
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Reports/Copies of SWS Census reports/Copies of SWS Census reports/2015/georeferenced"
transect2015<-readOGR(tdir,layer="TransectLines2015")
tdir<-"C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Reports/Copies of SWS Census reports/Copies of SWS Census reports/2014/georeferenced"
transect2014<-readOGR(tdir,layer="TransectLines2014")

#give them line numbers
transect2014$Id<-as.character(transect2014$Id)
transect2014$Id[1:16]<-16:1
transect2014$Id[17:24]<-17:24

transect2015$Id<-as.character(transect2014$Id)
transect2015$Id[1:16]<-16:1
transect2015$Id[17:24]<-17:24

transect2016$Id<-as.character(transect2016$Id)
transect2016$Id[1:16]<-16:1
transect2016$Id[17:24]<-17:24

#plot them
plot(myGrid3km)
plot(transect2014,col="blue",add=T)
plot(transect2015,col="red",add=T)
plot(transect2016,col="black",add=T)

# Intersect lines with raster "polygons" and add length to new lines segments
#"+proj=longlat +datum=WGS84"
#"+proj=eck4"
#2016
library(rgeos)
myGridproj<-myGrid3km
transect2016proj<-spTransform(transect2016,crs(myGrid3km))
rsp <- rasterToPolygons(myGridproj)
rp <- intersect(transect2016proj, rsp)
rp <- spTransform(rp,CRS("+proj=eck4"))
rp$length <- gLength(rp, byid=TRUE) 
rp2016<-rp
x <- tapply(rp$length, rp$layer, sum)
r2016 <- raster(myGrid3km)
r2016[as.integer(names(x))] <- x
plot(r2016)
#we now have the length that each transect traverses each grid cell by
#2015
transect2015proj<-spTransform(transect2015,crs(myGrid3km))
rsp <- rasterToPolygons(myGridproj)
rp <- intersect(transect2015proj, rsp)
rp <- spTransform(rp,CRS("+proj=eck4"))
rp$length <- gLength(rp, byid=TRUE) 
rp2015<-rp
x <- tapply(rp$length, rp$layer, sum)
r2015 <- raster(myGrid3km)
r2015[as.integer(names(x))] <- x
plot(r2015)

#2014
transect2014proj<-spTransform(transect2014,crs(myGrid3km))
rsp <- rasterToPolygons(myGridproj)
rp <- intersect(transect2014proj, rsp)
rp <- spTransform(rp,CRS("+proj=eck4"))
rp$length <- gLength(rp, byid=TRUE) 
rp2014<-rp
x <- tapply(rp$length, rp$layer, sum)
r2014 <- raster(myGrid3km)
r2014[as.integer(names(x))] <- x
plot(r2014)

#add to dataset

#there are 156 cells
length(getValues(r2014))==length(getValues(myGrid3km))
nrow(myGridDF3km)
myGridDF3km<-myGridDF3km[order(myGridDF3km$Grid3km),]
myGridDF3km$t2014<-getValues(r2014)
myGridDF3km$t2015<-getValues(r2015)
myGridDF3km$t2016<-getValues(r2016)
myGridDF3km$t2014[is.na(myGridDF3km$t2014)]<-0
myGridDF3km$t2015[is.na(myGridDF3km$t2015)]<-0
myGridDF3km$t2016[is.na(myGridDF3km$t2016)]<-0

#get information on cells in relationship to transect number
rp<-rbind(rp2016[,c("Id","layer")],
          rp2015[,c("Id","layer")],
          rp2014[,c("Id","layer")])
rp<-unique(data.frame(rp))
names(rp)[2]<-"Grid3km"

#add the digits together of the overlap transects
transectConversion<-ddply(rp,.(Grid3km),summarise,
                          transectId=paste(Id,collapse="."))

#make some changes to make better groupings
transectConversion$transectId[which(transectConversion$Grid3km==23)]<-"16.15"
transectConversion$transectId[which(transectConversion$Grid3km==31)]<-"15.14"
transectConversion$transectId[which(transectConversion$Grid3km==89)]<-"18.19.20"
transectConversion$transectId[which(transectConversion$Grid3km==95)]<-"4.3"
transectConversion$transectId[which(transectConversion$Grid3km==104)]<-"2.1"

myGridDF3km$transectId<-transectConversion$transectId[match(myGridDF3km$Grid3km,
                                                            transectConversion$Grid3km)]
myGridDF3km$transectId<-as.numeric(as.factor(myGridDF3km$transectId))
myGridDF3km$transectId[is.na(myGridDF3km$transectId)]<-max(myGridDF3km$transectId,na.rm=T)+1

#merge positive observation and transect effort

#call datafile
transectLengths<-data.frame(myGridDF3km[,c("Grid3km","t2014","t2015","t2016")])
transectLengths<-subset(transectLengths,rowSums(transectLengths[,2:4])>0)
transectLengths<-melt(transectLengths,id="Grid3km")  
transectLengths$Year<-as.numeric(gsub("t","",as.character(transectLengths$variable)))
names(transectLengths)[which(names(transectLengths)=="value")]<-"transectLength"

#merge effort with observations
transectDF$cells[!transectDF$cell %in% transectLengths$Grid3km]#none
datafile<-merge(transectLengths,transectDF,by.x=c("Grid3km","Year"),by.y=c("cells","Year"),all=T)
datafile$Count[is.na(datafile$Count)]<-0 

#only use data from our cells of interest
datafile<-subset(datafile,Grid3km%in%grid3kmImpute)
datafile<-subset(datafile,!is.na(transectLength))
nrow(datafile)

#convert grid to a continuous transect number
gridTranslate<-data.frame(Grid=sort(unique(datafile$Grid3km)))
gridTranslate$Transect<-1:nrow(gridTranslate)

#scale vars
datafile$Transect<-gridTranslate$Transect[match(datafile$Grid3km,gridTranslate$Grid)]
datafile$T<-datafile$Year-min(datafile$Year)+1
datafile$Obs<-ifelse(datafile$Count==0,0,1)

#re-ordering
datafile<-datafile[order(datafile$Year),]
datafile<-datafile[order(datafile$Transect),]
datafile<-datafile[order(datafile$Obs,decreasing=T),]
sites.lt<-gridTranslate$Grid  

#TransectYears
transectInfo<-unique(datafile[,c("Transect","T")])

#Get number of groups seen per transect
groupInfo<-acast(datafile,Transect~T,value.var="Obs",fun=sum,na.rm=T)

#Get total individuals seen per transects
totalsInfo<-acast(datafile,Transect~T,value.var="GroupSize",fun=sum,na.rm=T)

#get average group size per transect
groupSizes<-acast(datafile,Transect~T,value.var="GroupSize",fun=mean)

#transect lengths
mytransectLengths<-acast(datafile,Transect~T,value.var="transectLength",fun=mean)
mytransectLengths[is.na(mytransectLengths)]<-0

#detectionInfo
detectionInfo<-datafile[,c("Grid3km","Transect","Obs","Distance","GroupSize","T")]
detectionInfo<-subset(detectionInfo,Obs==1)
detectionInfo$Site[1]<-1
for(i in 2:nrow(detectionInfo)){
  detectionInfo$Site[i]<-ifelse(detectionInfo$Transect[i]==detectionInfo$Transect[i-1],detectionInfo$Site[i-1],
                                detectionInfo$Site[i-1]+1)
}

#detection covariates
detectionInfo$military<-myGridDF3km$Military[match(detectionInfo$Grid3km,myGridDF3km$Grid3km)]
detectionInfo$forestcover<-myGridDF3km$ForestCover[match(detectionInfo$Grid3km,myGridDF3km$Grid3km)]

#add transect grouping variable
detectionInfo$transectID<-myGridDF3km$transectId[match(detectionInfo$Grid3km,myGridDF3km$Grid3km)]
detectionInfo$transectID<-as.numeric(as.factor(detectionInfo$transectID))

#using RDistance
par(mfrow=c(1,1))
library(Rdistance)
fit <- F.dfunc.estim(detectionInfo$Distance[!is.na(detectionInfo$Distance)], likelihood="halfnorm")#sigma is 125.627 
plot(fit)

#index all data points
datafileObs<-subset(datafile,Obs==1)
TransYrIdx<-matrix(nrow=nrow(datafileObs),ncol=nrow(transectInfo))
TransYrIdx[]<-0
for(i in 1:nrow(datafileObs)){
  TransYrIdx[i,which(transectInfo$T==datafileObs$T[i]&
                       datafile$Transect[i]==transectInfo$Transect)]<-1
}

#################################################################################################

#get site indices
head(myGridDF3km)

#subset to those rows with data
myGridDF3km<-subset(myGridDF3km,Grid3km%in%c(sites.ct,sites.lt))

#relabel the site id
myGridDF3km$Grid<-1:nrow(myGridDF3km)
sites.lt<-sapply(sites.lt,function(x)myGridDF3km$Grid[which(myGridDF3km$Grid3km==x)])
sites.ct<-sapply(sites.ct,function(x)myGridDF3km$Grid[which(myGridDF3km$Grid3km==x)])

#################################################################################################

#investigate association between covariates
hist(myGridDF3km$ForestCover)
hist(myGridDF3km$Fields)
hist(myGridDF3km$Villages)
cor.test(myGridDF3km$Fields,myGridDF3km$ForestCover)#-0.3724837  
cor.test(myGridDF3km$Villages,myGridDF3km$Fields)#0.7962286  
cor.test(log(myGridDF3km$Villages+1),log(myGridDF3km$Fields+1))#0.8016456 
cor.test(myGridDF3km$Military,log(myGridDF3km$Villages+1))#0.2518769

##########################################################################################

#get centroids of the grid cells
#plot(myGrid3km)
tempDF<-as.data.frame(myGrid3km,xy=T)
myGridDF3km<-merge(myGridDF3km,tempDF,by.x="Grid3km",by.y="layer")

###########################################################################################
par(mfrow=c(1,1))
#plotting covariates
check<-rasterFromXYZ(myGridDF3km[,c("x","y","Military")])
plot(sws)
plot(check,col=gray.colors(n=10,start=0.9,end=0.3),add=T,
     axis.args=list(cex.axis=0.8),
     legend.args=list(text='Military area cover', side=4, line=2.5, cex=1.5))
plot(sws,add=T)

check<-rasterFromXYZ(myGridDF3km[,c("x","y","ForestCover")])
plot(sws)
plot(check,col=gray.colors(n=20,start=0.9,end=0.3),add=T,
    axis.args=list(cex.axis=0.8),
    legend.args=list(text='Forest cover', side=4, line=2.5, cex=1.5))
plot(sws,add=T)

check<-rasterFromXYZ(myGridDF3km[,c("x","y","Villages")])
plot(sws)
plot(check,col=gray.colors(n=20,start=0.9,end=0.3),add=T,
    axis.args=list(cex.axis=0.8),
    legend.args=list(text='Village pressure', side=4, line=2.75, cex=1.5))
plot(sws,add=T)

check<-rasterFromXYZ(myGridDF3km[,c("x","y","Fields")])
plot(sws)
plot(check,col=gray.colors(n=20,start=0.9,end=0.3),add=T,
    axis.args=list(cex.axis=0.8),
    legend.args=list(text='Field cover', side=4, line=2.5, cex=1.5))
plot(sws,add=T)

#################################################################################################
