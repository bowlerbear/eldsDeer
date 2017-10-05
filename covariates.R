#Looking at the patterns in spatial covariates across the sanctuary from remote sensing datasets

#libraries needed
library(rgdal)
library(raster)

################
#get study area#
################

#get shape file of the Sanctuary
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/ProtectedArea")
sws<-readOGR(dsn=getwd(),layer="WDPA_Aug2017_protected_area_1236-shapefile-polygons")

#look at the extent of the study region
sws@bbox

#get location of transects
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated")
transects<-readOGR(getwd(),layer="Transects")

#overlay transects with plot of the study area
plot(sws)
plot(transects,add=T)

#dummy raster
r<-raster(extent(sws@bbox))
res(r)<-0.01
r<-setValues(r,rnorm(ncell(r)))

########################
#get data on population#
########################

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers")
villages<-readOGR(dsn=getwd(),layer="mmr_mgy_pplp2_250k_wfp-mimu")

#overlay with larger empty raster
rL<-raster(extent(c(94.3,94.8,20,20.4)))
res(rL)<-0.01
rL<-setValues(rL,rnorm(ncell(rL)))
rasterCoords<-as.data.frame(xyFromCell(rL, 1:ncell(rL)))#get coordinates of the raster

#for each combination of centroids, get distances to all villages
out<-distanceFromPoints(rL,as(villages[,c("Longitude","Latitude")],'SpatialPoints'))
plot(out,main="Distance from villages")
plot(sws,add=T)

################
#get other data#
################

############
#topography#
############

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated")
top<-raster("MMR_alt.grd")
#res(top)
#[1] 0.008333333 0.008333333
#extract values for cells of our raster
topoSWS<-as.data.frame(extract(top,rasterCoords,cellnumbers=TRUE,na.rm=T))
rasterCoords$topo<-topoSWS$MMR_alt
topoRaster<-rasterize(rasterCoords[,1:2],r,field=rasterCoords[,3])
plot(topoRaster,main="Topography")
plot(sws,add=T)

############
#land cover#
############

landcover<-raster("MMR_cov.grd")
landcoverSWS<-as.data.frame(extract(landcover,rasterCoords[,1:2],cellnumbers=TRUE,na.rm=T))
rasterCoords$landcover<-landcoverSWS$MMR_cov
landcoverRaster<-rasterize(rasterCoords[,1:2],r,field=rasterCoords[,4])
plot(landcoverRaster,main="Basic Landcover")
plot(landcoverRaster,col=c("darkgreen","gray","yellow"))
plot(sws,add=T)
#from GLC2000
table(getValues(landcoverRaster))#most are 2 and 16

###########
#GlobCover#
###########

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated/Land cover, southern Asia (GlobCover 2009)")
globCover<-raster("layer11.tif")
sws2<-spTransform(sws,CRS(projection(globCover)))
globCover<-crop(globCover,extent(sws2@bbox))
library(RColorBrewer)
plot(globCover,col=sample(brewer.pal(11,"Set3")))
plot(sws2,add=T)

#############
#Hansen data#
############

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated/Hansen")
forest2000<-raster("Hansen_GFC2015_treecover2000_30N_090E.tif")
forest2000<-crop(forest2000,extent(sws@bbox))
plot(forest2000,main="ForestCover in 2000")
plot(sws,add=T)

##############
#Chelsea data#
##############

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated/CHELSA_tmax10_3_land")
temp<-raster("CHELSA_tmax10_3_1979-2013_V1.2_land.tif")
temp<-crop(temp,extent(sws@bbox))
plot(temp,main="max temp in March")
plot(sws,add=T)

#oTHERS

#MODIS data
#climatic data - use CRU?
#temporal data or spatial
