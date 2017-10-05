
#Get the formatted data for BUGS
source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/GitFiles/formattingCameratraps_occupancyModel.R')

#####################################################################################

#basic site occupancy including detection and occupancy covariates

str( win.data <- list(y = y, M = nrow(y), J = ncol(y), siteCovs = occupancyCovariates[,2:3],
                      WATER = WATER, YEAR = YEAR) )

# Specify model in BUGS language
sink("model.txt")
cat("
    model {
    
    # Priors
    mean.psi ~ dunif(0, 1)
    mean.p ~ dunif(0, 1)

    #priors for occupancy covariates
    beta0 <- logit(mean.psi)
    beta.forest ~ dnorm(0,0.0001)
    beta2.forest ~ dnorm(0,0.0001)

    #priors for detection 
    alpha0 <- logit(mean.p)
    alpha.water ~ dnorm(0,0.0001)
    alpha.year ~ dnorm(0,0.0001)


    # Likelihood
    for (i in 1:M) {    # Loop over sites
    # State model
      z[i] ~ dbern(psi[i])         
      logit(psi[i]) <- beta0 + beta.forest*siteCovs[i,1] + beta2.forest*siteCovs[i,2]
    }

    for (i in 1:M) {
      for (j in 1:J) { # Loop over replicate surveys
      y[i,j] ~ dbern(z[i]*p[i,j])  
      # Observation model
      logit(p[i,j]) <- alpha0 + alpha.water*WATER[i,j] + alpha.year*YEAR[i,j]
      #logit(p[i,j]) <- alpha0
    }
    }

    }
    ",fill = TRUE)
sink()

# Initial values
zst <- apply(y, 1, max, na.rm=T)  
inits <- function(){list(z = zst, mean.psi = rnorm(1,0.3,0.05), mean.p = rnorm(1,0.2,0.05))}

# Parameters monitored
params <- c("mean.p","alpha.water","alpha.year","mean.psi","beta.forest","beta2.forest","psi","z","p")

# MCMC settings
ni <- 10000   ;   nt <- 1   ;   nb <- 1000   ;   nc <- 3

# Call JAGS and summarize posteriors
library(jagsUI)
fm2 <- jags(win.data, inits, params, "model.txt", n.chains = nc,
   n.thin = nt, n.iter = ni, n.burnin = nb)
print(fm2, dig = 3)

#compare, z, psi and zst
out<-data.frame(Psi=fm2$mean$psi,Z=fm2$mean$z,zst)

#with constant p
summary(subset(out,zst==1))
Psi                Z          zst   
Min.   :0.05271   Min.   :1   Min.   :1  
1st Qu.:0.37992   1st Qu.:1   1st Qu.:1  
Median :0.37992   Median :1   Median :1  
Mean   :0.35714   Mean   :1   Mean   :1  
3rd Qu.:0.37992   3rd Qu.:1   3rd Qu.:1  
Max.   :0.39425   Max.   :1   Max.   :1

summary(subset(out,zst==0))
Psi                 Z                  zst   
Min.   :0.000942   Min.   :0.0000000   Min.   :0  
1st Qu.:0.002189   1st Qu.:0.0001111   1st Qu.:0  
Median :0.064721   Median :0.0024815   Median :0  
Mean   :0.169577   Mean   :0.0299151   Mean   :0  
3rd Qu.:0.379920   3rd Qu.:0.0397778   3rd Qu.:0  
Max.   :0.395716   Max.   :0.1906296   Max.   :0

###################################################################################

#Plotting the model predictions

#extract predicted psi for each grid cell
occupancyCovariates$PredictedPSI<-fm2$mean$psi

#get coordinates of each grid cell
gridCoords<-unique(locations[,c("GRID.CELL.ID","LATITUDE..N.","LONGITUDE..E.")])
gridCoords<-subset(gridCoords,!duplicated(GRID.CELL.ID))
gridCoords$PredictedPSI<-occupancyCovariates$PredictedPSI[match(gridCoords$GRID.CELL.ID,occupancyCovariates$GRID.CELL.ID)]

#plot the data
coordinates(gridCoords)<-c("LONGITUDE..E.","LATITUDE..N.")
overlapCell<-data.frame(extract(r,gridCoords,cellnumbers=T))
gridCoords$cell<-as.numeric(overlapCell$cells)
rasterDF<-data.frame(xyFromCell(r,1:ncell(r)))
rasterDF$cellnu<-as.numeric(cellFromXY(r,rasterDF))
rasterDF$PredictedPSI<-gridCoords$PredictedPSI[match(rasterDF$cellnu,gridCoords$cell)]
predictedRaster<-rasterFromXYZ(rasterDF[,c("x","y","PredictedPSI")])
plot(predictedRaster)

#predict to whole region on the basis of forect cover
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/MMR layers/Updated/Hansen")
forest2000<-raster("Hansen_GFC2015_treecover2000_30N_090E.tif")
forest2000<-crop(forest2000,extent(r))
#project this raster onto the grid cell of the camera traps
forest2000DF<-as.data.frame(forest2000,xy=T)
coordinates(forest2000DF)<-c("x","y")
forestOverlap<-data.frame(extract(r,forest2000DF,cellnumbers=T))
forest2000DF$cell<-as.numeric(forestOverlap$cells)
forestOverlap<-ddply(forest2000DF@data,.(cell),summarise,ForestCover=mean(Hansen_GFC2015_treecover2000_30N_090E,na.rm=T))
rasterDF$ForestCover<-forestOverlap$ForestCover[match(rasterDF$cellnu,forestOverlap$cell)]

#make predictions from the regression model
library(boot)
rasterDF$MyPredictions<-inv.logit(logit(fm2$mean$mean.psi)+fm2$mean$beta.forest*rasterDF$ForestCover+fm2$mean$beta2.forest*rasterDF$ForestCover^2)

#plot them
predictedRaster <- crop(predictedRaster, extent(sws))
predictedRaster <- mask(predictedRaster, sws)
plot(predictedRaster)
plot(sws,add=T)

######################################################################