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
locations$Site<-paste(locations$LONGITUDE..E.,locations$LATITUDE..N.)
surveysRA<-subset(locations,!Site%in%surveysP@data$Site)
coordinates(surveysRA)<-c("LONGITUDE..E.","LATITUDE..N.")
#identify the cells
r[]<-1:ncell(r)
surveysRA<-extract(r,surveysRA)

#from the transects
transectsRA<-do.call(c,extract(r,transects))
transectsRA<-transectsRA[!transectsRA%in%transectPointsRP]

#plotting presences and absences
r[]<-0
r[surveysRA]<-1
r[transectsRA]<-2
r[surveysRP]<-4
r[transectPointsRP]<-5
plot(r,col=c("white","azure2","azure3","pink","red"))
plot(sws,add=T)

#######################
#Obtain the covariates#
#######################

#is a cell inside or outside the military area
r[]<-1:ncell(r)
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers")
military<-readOGR(getwd(),layer="Military area")
militaryR<-extract(r,military)[[1]]
#plot(sws)
#plot(military,add=T,col="grey")

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

#import Johns work in google earth
fields<-readOGR(getwd(),layer="JohnsPolygons")
fields<-subset(fields,grepl("Field",fields$Name))
fields<-subset(fields,!grepl("military",fields$Name))
#plot(sws)
#plot(fields,add=T,col="grey")
#calculate coverage of each grid cell by a field
fieldsR<-ldply(extract(r,fields,weights=T,normalizeWeights=F))
fieldsR<-ddply(fieldsR,.(value),summarise,weight=sum(weight))

#also get information on nearby villages
villages<-readOGR(getwd(),layer="JohnsPolygons")
villages<-subset(villages,grepl("Village",villages$Name))
plot(sws)
plot(villages,add=T,col="grey")
#extract density
#probably should increase the extent of r
r2<-r
extent(r2)<-c(94.2,94.9,20,20.4)
res(r2)<-res(r)
r2[]<-1:ncell(r2)
villagesR<-ldply(extract(r2,villages,weights=T,normalizeWeights=F))
#sum it up over each cell
villagesR<-ddply(villagesR,.(value),summarise,weight=sum(weight))
#overlay onto the raster
rVillages<-as.data.frame(r2,xy=T)
rVillages$Village<-villagesR$weight[match(rVillages$layer,villagesR$value)]
rVillages$Village[is.na(rVillages$Village)]<-0
villageRaster<-rasterFromXYZ(rVillages[,c("x","y","Village")])

library(gdistance)

#smooth the raster
villageRaster <- focal(villageRaster, w=matrix(1, 17, 17), mean,na.rm=T)
projection(villageRaster)<-CRS("+proj=longlat +ellps=WGS84") 
villageRaster<-projectRaster(from=villageRaster,to=myGrid)
rVillages<-as.data.frame(r2,xy=T)
plot(villageRaster)
plot(sws,add=T)
plot(villages,add=T)
villagesDF<-as.data.frame(villageRaster)
villagesDF$cell<-1:ncell(villageRaster)

##################
#combine all data#
##################

#on the 1km grid
myGridDF$Military<-sapply(myGridDF$Grid1km,function(x)ifelse(x%in%militaryR,1,0))
myGridDF$ForestCover<-forestDF$ForestCover[match(myGridDF$Grid1km,forestDF$Cell)]
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
                   Military=sum(Military),
                   Villages=sum(Villages,na.rm=T),
                   ForestCover=sum(ForestCover),
                   Fields=sum(Fields),
                   Deer.ct=max(Deer.ct,na.rm=T),
                   Deer.lt=max(Deer.lt,na.rm=T),
                   TransectNu=sum(Transect,na.rm=T))
myGridDF3km$Deer.ct[is.infinite(myGridDF3km$Deer.ct)]<-NA
myGridDF3km$Deer.lt[is.infinite(myGridDF3km$Deer.lt)]<-NA
myGridDF3km$Deer<-ifelse(myGridDF3km$Deer.ct==1|myGridDF3km$Deer.lt==1,1,0)

#also add some lat and lon coordinates

########################################################################

#Fitting a hierarcical model to both datasets

#get full list of grid cells with data

##########################
#(1) get camera trap data#
##########################

#initally use the 1km grid
cameraTraps<-merge(surveys,locations,by="ID")
coordinates(cameraTraps)<-c("LONGITUDE..E.","LATITUDE..N.") 
cameraTraps$Cell<-extract(myGrid,cameraTraps)
cameraTraps<-merge(myGridDF,data.frame(cameraTraps),by.x="Grid1km",by.y="Cell")
cameraTraps$Year<-sapply(as.character(cameraTraps$ID),function(x)strsplit(x,"_")[[1]][1])

#get rid of data without water information
nowaterData<-locations$ID[locations$WATER.POINT==""]
cameraTraps<-subset(cameraTraps,!ID%in%nowaterData)

#combine water point or salt lick as an attractor
attractorDF<-unique(cameraTraps[,c("ID","Grid1km","WATER.POINT","SALT.LICK")])

#check attractor is constant within a grid cell
out<-ddply(attractorDF,.(Grid1km),summarise,nu=length(unique(WATER.POINT)))
#5 still have variation in whether it is a water point
dups<-subset(attractorDF,Grid1km%in%out$Grid1km[out$nu==2])

#first drop those with 'no' water points
cameraTraps<-subset(cameraTraps,!cameraTraps$ID%in%dups$ID[dups$WATER.POINT=="No"])

#recount Day so it goes over multiple years/cameras in the same grid cells
cameraTraps<-ddply(cameraTraps,.(Grid1km),function(x){
  Day=1:nrow(x)
  data.frame(cbind(Day2=Day,x))
})

#relabel Grid to be 1 to 9
cameraTraps<-ddply(cameraTraps,.(Grid3km),function(x){
  x$GridRep=as.factor(x$Grid1km)
  levels(x$GridRep)<-1:length(unique(x$GridRep))
  x$GridRep<-as.numeric(as.character(x$GridRep))
  return(x)
})

#convert data into a 3-dimensional array
obsMatrix<-acast(cameraTraps,Grid3km~GridRep~Day2,value.var="PA")
dim(obsMatrix)

nu<-30
y<-obsMatrix[,,1:nu]

#water info
cameraTraps$WATER.POINT<-sapply(cameraTraps$WATER.POINT,function(x)ifelse(x=="Yes",1,0))
waterMatrix<-acast(cameraTraps,Grid3km~GridRep,value.var="WATER.POINT",fun=max)
waterMatrix[is.infinite(waterMatrix)]<-0

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
plot(transects,add=T)
plot(transectPoints,add=T,col="red")

#on point does not overlap on transect 6  -move it to the next grid 81
transectPoints$cells[transectPoints$cells==80]<-81

#get list of all cells that each transect overlapped
transectGrids<-do.call(c,extract(myGrid3km,transects))
transectGrids<-data.frame(Grid=transectGrids)
#expand to include sampling for each year
transectGrids<-expand.grid(Year=unique(transectPoints$Year),Grid=unique(transectGrids$Grid))

#merge
transectDF<-merge(as.data.frame(transectGrids),transectPoints@data,all=T,by.x=c("Grid","Year"),by.y=c("cells","Year"))
#Count is 1 when something was seen, 0 if not
transectDF$GroupSize<-transectDF$Count
transectDF$Obs[is.na(transectDF$Count)]<-0
transectDF$Obs[transectDF$Count!=0]<-1

#truncate observation at 250 m
transectDF$Count[transectDF$Distance>250]<-0
transectDF$Obs[transectDF$Distance>250]<-0
transectDF$GroupSize[transectDF$Distance>250]<-NA
transectDF$GroupSize[transectDF$GroupSize>250]<-NA

head(transectDF)

#call datafile
datafile<-transectDF

#convert grid to a continuous transect number
gridTranslate<-data.frame(Grid=sort(unique(datafile$Grid)))
gridTranslate$Transect<-1:nrow(gridTranslate)

#scale vars
datafile$Transect<-gridTranslate$Transect[match(datafile$Grid,gridTranslate$Grid)]
datafile$T<-datafile$Year-min(datafile$Year)+1

#get info on transect lengths for each year and transect
#transectDistances<-acast(datafile,Transect~T,value.var="TransectDist",fun=max)

#TransectYears
transectInfo<-unique(datafile[,c("Transect","T")])

#re-ordering
datafile<-datafile[order(datafile$Year),]
datafile<-datafile[order(datafile$Transect),]
datafile<-datafile[order(datafile$Obs,decreasing=T),]

#to align covariates
covAlign<-acast(datafile,Grid~T,value.var="Obs",fun=sum,na.rm=T)
sites.lt<-as.numeric(dimnames(covAlign)[[1]])  

#Get number of groups seen per transect
groupInfo<-acast(datafile,Transect~T,value.var="Obs",fun=sum,na.rm=T)

#Get total individuals seen per transects
totalsInfo<-acast(datafile,Transect~T,value.var="GroupSize",fun=sum,na.rm=T)

#get average group size per transect
groupSizes<-acast(datafile,Transect~T,value.var="GroupSize",fun=mean)

#detectionInfo
detectionInfo<-datafile[,c("Grid","Transect","Obs","Distance","GroupSize")]
detectionInfo<-subset(detectionInfo,Obs==1)
detectionInfo$Site[1]<-1
for(i in 2:nrow(detectionInfo)){
  detectionInfo$Site[i]<-ifelse(detectionInfo$Transect[i]==detectionInfo$Transect[i-1],detectionInfo$Site[i-1],
                                detectionInfo$Site[i-1]+1)
}


#detection covariates
detectionInfo$military<-myGridDF3km$Military[match(detectionInfo$Grid,myGridDF3km$Grid3km)]
detectionInfo$forestcover<-myGridDF3km$ForestCover[match(detectionInfo$Grid,myGridDF3km$Grid3km)]

#look at frequency histogram of distribution distances
library(lattice)
histogram(detectionInfo$Distance)

#add half normal line
x <- detectionInfo$Distance
origx<-x
x<-c(x,x*-1)
h<-hist(x, breaks=16) 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=0,sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2) 

h<-hist(origx, breaks=8, col="lightgray",ylim=c(0,55),xlab="Distance(m)") 
lines(xfit[xfit>0], yfit[xfit>0], col="blue", lwd=2)

#using RDistance
library(Rdistance)
fit <- F.dfunc.estim(origx, likelihood="halfnorm")#sigma is 125.627 
plot(fit)

#index all data points
datafileObs<-subset(datafile,Obs==1)
TransYrIdx<-matrix(nrow=nrow(datafileObs),ncol=nrow(transectInfo))
TransYrIdx[]<-0
for(i in 1:nrow(datafileObs)){
  TransYrIdx[i,which(transectInfo$T==datafileObs$T[i]&
                       datafile$Transect[i]==transectInfo$Transect)]<-1
}

#get occupancy data
pa<-groupInfo
pa[pa>0]<-1
pa<-apply(pa,1,function(x)ifelse(all(x==1),1,0))


####################################################
#Organise covariates into a dataframe for the model#
####################################################

head(myGridDF3km)

#Camera data is first and then the line transects
myGridDF3km<-subset(myGridDF3km,!is.na(Deer.ct)|!is.na(Deer.lt))
forestcover<-myGridDF3km$ForestCover
forestcoverB<-ifelse(forestcover<40,0,1)
forestcover2<-log(forestcover+1)
military<-myGridDF3km$Military
fields<-myGridDF3km$Fields
fields<-sqrt(fields+1)
villages<-myGridDF3km$Villages
presenceRecord<-ifelse(myGridDF3km$Deer==1&!is.na(myGridDF3km$Deer),1,0)

#relabel the site id
myGridDF3km$Grid<-1:nrow(myGridDF3km)
sites.lt<-sapply(sites.lt,function(x)myGridDF3km$Grid[which(myGridDF3km$Grid3km==x)])
sites.ct<-sapply(sites.ct,function(x)myGridDF3km$Grid[which(myGridDF3km$Grid3km==x)])

#get transect effort

#using a buffer
library(rgeos)
transectsProj<-spTransform(transects,CRS("+proj=eck4"))
transectsProj_Buffered<-gBuffer(transectsProj,width=250)
transects_Buffered<-spTransform(transectsProj_Buffered,crs(myGrid3km))
plot(myGrid3km)
plot(transects_Buffered,add=T)
#get area of overlap
transectOverlap<-data.frame(extract(myGrid3km,transects_Buffered,weights=T,normalizeWeights=F))

#add to dataset
myGridDF$Overlap<-transectOverlap$weight[match(myGridDF$Grid3km,transectOverlap$value)]
overlapRaster<-rasterFromXYZ(myGridDF[,c("x","y","Overlap")])
plot(overlapRaster)
#scale overlap between 0 and 3...?
myGridDFOverlap<-subset(myGridDF,!is.na(Overlap))
myGridDFOverlap$Overlap<-myGridDFOverlap$Overlap*9
transectAreas<-myGridDFOverlap$Overlap[match(gridTranslate$Grid,myGridDFOverlap$Grid3km)]

#or based on number of 1km grids being sampled
transectAreas<-myGridDF3km$TransectNu[match(sites.lt,myGridDF3km$Grid)]
transectAreas[transectAreas==0]<-1
transectAreas<-transectAreas/9

#################################################################################################

#investigate association between covariates
cor.test(myGridDF3km$Villages,myGridDF3km$Fields)#0.4470563 
hist(sqrt(myGridDF3km$Villages+1))
hist(sqrt(myGridDF3km$Fields+1))
cor.test(sqrt(myGridDF3km$Villages+1),sqrt(myGridDF3km$Fields+1))#0.5031411 
#################################################################################################