#get the file that already has been partly formatted
source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/GitFiles/formattingCameratraps.R')

#needs fixing

######################
###multi-scale model##
#####################

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/Camera traps")

#get camera trap grid
grid<-read.delim("CameraGrid.txt")
grid<-subset(grid,Col!=19)
library(raster)
gridExtent<-c(94.35,94.69,19.98,20.36)
r<-raster(extent(gridExtent),nrow=21,ncol=18)#divisable by 3
r[]<-grid$Grid
projection(r)<-CRS("+proj=longlat +datum=WGS84")
plot(r)
rDF<-as.data.frame(r,xy=T)
rDF$cellnu<-as.numeric(cellFromXY(r,rDF[,1:2]))

#group cells into 9 by 9 grid
r2<-raster(extent(gridExtent),ncol=6,nrow=7)
r2[]<-1:ncell(r2)
#disaggregate r2
r2<-disaggregate(r2,fact=3)
r2DF<-as.data.frame(r2,xy=T)
r2DF$cellnu<-as.numeric(cellFromXY(r2,r2DF[,1:2]))
rDF$LargeCell<-r2DF$layer

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

#add 9 cell grouping to the surveys
surveys$SimpleGrid<-gsub("SWS","",surveys$Grid)
surveys$SimpleGrid<-as.numeric(surveys$SimpleGrid)
surveys$SimpleGrid<-rDF$cellnu[match(surveys$SimpleGrid,rDF$layer)]
surveys$SimpleGrid[which(surveys$Grid=="RMY")]<-166#the army grid cell is cell number 166
surveys$LargeCell<-rDF$LargeCell[match(surveys$SimpleGrid,rDF$cellnu)]

#relabel Grid to be 1 to 9
surveys<-ddply(surveys,.(LargeCell),function(x){
  x$GridRep=as.factor(x$SimpleGrid)
  levels(x$GridRep)<-1:length(unique(x$GridRep))
  x$GridRep<-as.numeric(as.character(x$GridRep))
  return(x)
})


#convert data into a 3-dimensional array
obsMatrix<-acast(surveys,LargeCell~GridRep~Day2,value.var="PA")
dim(obsMatrix)#21,9 and 99

nu<-30
y<-obsMatrix[,,1:nu]

#9-unit square
#square

#occupancy factors
Covariates<-unique(locations[,c("GRID.CELL.ID","ForestCover")])
Covariates$SimpleGrid<-as.numeric(gsub("SWS","",Covariates$GRID.CELL.ID))#army cell gets NA for the moment
Covariates$SimpleGrid<-rDF$cellnu[match(Covariates$SimpleGrid,rDF$layer)]
Covariates$SimpleGrid[which(Covariates$Grid=="RMY")]<-166
Covariates$LargeCell<-rDF$LargeCell[match(Covariates$SimpleGrid,rDF$cellnu)]

occupancyCovariates<-ddply(Covariates,.(LargeCell),summarise,ForestCover=mean(ForestCover))
occupancyCovariates$ForestCoverSQD<-(occupancyCovariates$ForestCover+0.001)^2
occupancyCovariates[,2:3]<-sapply(occupancyCovariates[,2:3],scale)

#use covariates
useCovariates<-ddply(Covariates,.(LargeCell,SimpleGrid),summarise,ForestCover=mean(ForestCover))
useCovariates$ForestCoverSQD<-(useCovariates$ForestCover+0.001)^2
useCovariates[,3:4]<-sapply(useCovariates[,3:4],scale)
useCovariates$GridRep<-surveys$GridRep[match(useCovariates$SimpleGrid,surveys$SimpleGrid)]
useMatrix_Forest<-acast(useCovariates,LargeCell~GridRep,value.var="ForestCover")
useMatrix_ForestSQD<-acast(useCovariates,LargeCell~GridRep,value.var="ForestCoverSQD")

#detection covariates
detectionCovariates<-unique(locations[,c("ID","GRID.CELL.ID","Year","WATER.POINT")])
detectionCovariates<-merge(surveys,detectionCovariates,by=c("ID","Year"),all.x=T)
detectionCovariates$Year<-as.numeric(detectionCovariates$Year)
detectionCovariates$WATER.POINT<-ifelse(detectionCovariates$WATER.POINT=="Yes",1,0)
WATER<-acast(detectionCovariates,LargeCell~GridRep~Day2,value.var="WATER.POINT")[,,1:nu]
YEAR<-acast(detectionCovariates,LargeCell~GridRep~Day2,value.var="Year")[,,1:nu]

#put any value for the missing data
WATER[is.na(WATER)]<-0
YEAR[is.na(YEAR)]<-2016
YEAR[YEAR==2014]<-0
YEAR[YEAR==2016]<-1

dim(WATER)
dim(YEAR)

############################################################################
#basic site occupancy including detection and occupancy covariates

#with 36 km2
str( win.data <- list(y = y, n.largecell = dim(y)[1], n.smallcell = dim(y)[2], n.rep = dim(y)[3], 
                      largecellCovs = occupancyCovariates[2:3], 
                      useF = useMatrix_Forest, useF2 = useMatrix_ForestSQD, 
                      WATER = WATER, YEAR = YEAR))

#with 4 km2 cell
y<-obsMatrix
str( win.data <- list(y = y, n.largecell = dim(y)[1], n.smallcell = dim(y)[2], n.rep = dim(y)[3], 
                      largecellCovs = occupancyCovariates, 
                      WATER = WATER, YEAR = YEAR))


# Specify model in BUGS language
sink("model.txt")
cat("
    model {
    
    # Priors and model for params
    
    # Intercept of occupancy probability
    int.psi ~ dunif(0,1)         
    
    # Intercepts availability probability
    for(t in 1:n.smallcell){
    int.theta[t] ~ dunif(0,1) 
    }

    # Intercepts detection probability
    for(t in 1:n.rep){
    int.p[t] ~ dunif(0,1)     
    }
    
    #covariates
    occ.forest ~ dnorm(0, 0.1)
    occ.forest2 ~ dnorm(0, 0.1)    
    #avail.forest ~ dnorm(0, 0.1)
    #avail.forest2 ~ dnorm(0, 0.1)
    avail.water ~ dnorm(0, 0.1)
    det.year ~ dnorm(0, 0.1)
    
    # 'Likelihood' 
    for (i in 1:n.largecell){
    z[i] ~ dbern(psi[i])
    logit(psi[i]) <- logit(int.psi) + occ.forest*largecellCovs[i,1] + occ.forest2*largecellCovs[i,2]
    
    for (j in 1:n.smallcell){
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z[i] * theta[i,j]
    logit(theta[i,j]) <- logit(int.theta[j]) 
    
    for (k in 1:n.rep){
    y[i,j,k] ~ dbern(mu.y[i,j,k])
    mu.y[i,j,k] <- a[i,j] * p[i,j,k]
    logit(p[i,j,k]) <- logit(int.p[k]) +  det.year * YEAR[i,j,k]+ avail.water * WATER[i,j,k]
    }
    }
    }
    
    } # end model
    ",fill=TRUE)
    sink()
    
# Initial values
zst <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

#inits <- function(){list(z = zst)}
inits <- function(){list(z = zst, a = ast,int.psi = rnorm(1,0.3,0.05))}

# Parameters monitored
params<-c("int.psi","occ.forest","occ.forest2","avail.water","det.year","psi")

# MCMC settings
ni <- 10000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3

# Call JAGS and summarize posteriors
library(jagsUI)
fm2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm2, dig = 3)

###################################################################################
