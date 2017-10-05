#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/GitFiles/formattingcombinedModel.R')

######################################################################################

#Fitting a hierarcical model to both datasets

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

head(myGridDF3km)

#Others
forestcover<-c(myGridDF3km$ForestCover[as.numeric(dimnames(y)[[1]])])
forestcover<-ifelse(forestcover<10,0,1)
fields<-c(myGridDF3km$Fields[as.numeric(dimnames(y)[[1]])])
fields<-sqrt(fields+1)
#villages<-c(myGridDF3km$Villages[as.numeric(dimnames(y)[[1]])])

########################
#Compile data for model#
########################

bugs.data<-list(#camera trap data
                y.ct = y,
                n.CameraTrapSites = dim(y)[1],
                n.reps = dim(y)[3],
                n.1kmGrids = dim(y)[2],
                #covariates
                Forest=forestcover,
                Field=fields,
                Villages=villages,
                Water=waterMatrix,
                Month=monthMatrix)


sink("combinedModelv2.txt")
cat("
    model {
    
    # Priors
    mean.psi ~ dunif(0, 1)
    mean.p.ct ~ dunif(0, 1)

    # Intercepts availability probability
    for(j in 1:n.1kmGrids){
    int.theta[j] ~ dunif(0,1) 
    }
   
    #priors for detection 
    alpha.ct <- logit(mean.p.ct)
    theta.water ~ dnorm(0,0.001)
    alpha.month ~ dnorm(0,0.001)
 
    #priors for occupancy 
    beta0 <- logit(mean.psi)
    beta.forest ~ dnorm(0,0.001)
    #beta.villages ~ dnorm(0,0.001)
    #beta.fields ~ dnorm(0,0.001)

    
    #(1) Camera trap data:
    for (i in 1:n.CameraTrapSites) { 
    logit(psi[i]) <- beta0 + beta.forest * Forest[i]
    z.ct[i] ~ dbern(psi[i])
  
    #availbility model
    for (j in 1:n.1kmGrids){
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z.ct[i] * theta[i,j]
    logit(theta[i,j]) <- logit(int.theta[j]) + theta.water * Water[i,j]
  
    for (k in 1:n.reps) { # Loop over replicate surveys
    y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
    mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
    logit(p.ct[i,j,k]) <- alpha.ct + alpha.month * Month[i,j,k]
    }
    }
    }
    
    }
    ",fill = TRUE)
sink()

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')


zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0


params <- c("beta0","beta.forest","theta.water","alpha.month","mean.p.ct","psi")

inits <- function(){list(z.ct = zst.ct,
                         a = ast)}

#n.iter<-10000

out1 <- jags(bugs.data, inits=inits, params, "combinedModelv2.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
traceplot(out1)

########################################################################################

#get occurence data

#just over the range od data
PredictedPSI<-out1$mean$psi
cameraTrapSites<-data.frame(Sites=as.numeric(dimnames(y)[[1]]))
cameraTrapSites$Preds<-PredictedPSI[1:bugs.data$n.CameraTrapSites]
myGridDF$Preds<-cameraTrapSites$Preds[match(myGridDF$Grid3km,cameraTrapSites$Sites)]
predRaster<-rasterFromXYZ(myGridDF[,c("x","y","Preds")])
plot(predRaster)

#imputing over the whole range
library(boot)
myGridDF<-subset(myGridDF,!is.na(ForestCover))
forestcover<-ifelse(myGridDF$ForestCover<10,0,1)
myGridDF$Preds<-inv.logit(out1$mean$beta0+out1$mean$beta.forest*forestcover)
predRaster<-rasterFromXYZ(myGridDF[,c("x","y","Preds")])
plot(predRaster)
predRaster<-mask(predRaster,sws)
plot(predRaster)
plot(sws,add=T)

########################################################################################
