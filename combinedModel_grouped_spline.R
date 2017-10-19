#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

######################################################################################

#get grid with simple fits on from the script of 'combinedModel_grouped'
load("myGridDF3km_wFits.RData")

#get centroids of the grid cells for spline
plot(myGrid3km)
tempDF<-as.data.frame(myGrid3km,xy=T)
myGridDF3km<-merge(myGridDF3km,tempDF,by.x="Grid3km",by.y="layer")

#Use the JAGAM function to get the BUGS code for the spline
#http://www.petrkeil.com/?p=2385
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
library(mgcv)
jags.ready <- jagam(round(fits)~s(x, y), 
                    data=myGridDF3km, 
                    family="poisson", 
                    file="jagam.txt")


#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero
  
#######################
#Compile data for model#
########################

bugs.data<-list(#camera trap data
                site.ct = sites.ct,
                y.ct = y,
                n.CameraTrapSites = dim(y)[1],
                n.reps = dim(y)[3],
                n.1kmGrids = dim(y)[2],
                #line transect data
                site.lt = sites.lt,
                n.Transect = length(unique(datafile$Transect)), 
                n.Yrs = length(unique(datafile$T)),
                W = 250,
                n = totalsInfo,
                n.Detections = nrow(detectionInfo),
                #n.detectionSites = length(unique(detectionInfo$Site)),
                #d.Forest = detectionInfo$forestcover,
                #d.Military = detectionInfo$military,
                y = detectionInfo$Distance,
                #d.Groupsize = detectionInfo$GroupSize,
                #d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectAreas = transectAreas,
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                #covariates
                #Forest=forestcoverB,
                #Forest2=forestcover2,
                #Fields=as.numeric(scale(sqrt(fields+1))),
                #Military=ifelse(military>0,1,0),
                #Villages=as.numeric(scale(sqrt(villages+1))),
                Water=waterMatrix,
                Month=monthMatrix,
                #spatial covariates
                X = X,
                S1 = S1,
                zero = zero)

#just using line transect data

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
sink("combinedModel_grouped_spline.txt")
cat("
    model {
    
    #Model of factors affecting abundance
    
    #the linear predictor
    eta <- X %*% b ## linear predictor
    
    for (i in 1:n.sites) { #across all sites   
      abund.effect[i] <-  eta[i]
    }

    #Model for the spline

    ## Parametric effect priors 
    ##CHECK tau=1/2^2 is appropriate!
    for (i in 1:1) { 
      #b[i] ~ dnorm(0,0.24) 
      b[i] ~ dnorm(0,0.001) 
    }
    
    ## prior for s(x,y)... 
    K1 <- S1[1:29,1:29] * lambda[1]  + S1[1:29,30:58] * lambda[2]
    b[2:30] ~ dmnorm(zero[2:30],K1) 
    
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
      #lambda[i] ~ dgamma(.05,.005)
        lambda[i] ~ dunif(0,1)
      rho[i] <- log(lambda[i])
    }

    # MODEL NUMBER OF INDIVIDUALS DETECTED

    # PRIORS
    intercept.lt ~ dunif(-10, 10)
    
    #fit model
    for (j in 1:n.Transect){
      for (t in 1:n.Yrs){
        n[j,t] ~ dpois(nHat[j,t]*transectAreas[j]*averagePa) 
        log(nHat[j,t]) <- abund.effect[site.lt[j]]
    }
    }

    #predict density (number of indivdiuals per grid cell) for all sites 
    #using the abundance effect and average group size
    #and the intercept of the line transect model for number of groups
    for(i in 1:n.sites){
      log(Density[i]) <- abund.effect[i]
    }

    # DETECTABILITY COMPONENT#####

    pi <- 3.141593
    
    #Priors
    #mean.p.lt ~ dunif(0, 1)
    #alpha.lt <- logit(mean.p.lt)
    sigma ~ dunif(50,140)
    b.d.0 ~ dunif(0,20)
    #b.d.GroupSize ~ dnorm(0,0.0001)#no effect
    #b.d.Forest ~ dnorm(0,0.0001)#no effect
    
    ##### Begin model for *all detections*
    
    for( i in 1:n.Detections){
    
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    mu.df[i] <- b.d.0
    
    # estimate of sd and var, given coeffs above
    sig.df[i] <- exp(mu.df[i])
    sig2.df[i] <- sig.df[i]*sig.df[i]
    
    # effective strip width
    esw[i] <- sqrt(pi * sig2.df[i] / 2)
    f0[i] <- 1/esw[i] #assuming all detected on the line
    
    # LIKELIHOOD
    # estimate sigma using zeros trick
    L.f0[i] <- exp(-y[i]*y[i] / (2*sig2.df[i])) * 1/esw[i] #y are the distances
    nlogL.f0[i] <-  -log(L.f0[i])
    zeros.dist[i] ~ dpois(nlogL.f0[i])
    }
    
    #get average esw per transect and year
    for(k in 1:n.TransectYrs){
    for(i in 1:n.Detections){
    grp.ESW[i,k] <- esw[i] * TransYrIdx[i,k]
    }
    ESW.jt[ty.combos[k,1], ty.combos[k,2]] <- sum(grp.ESW[,k])/max(1,sum(TransYrIdx[,k]))  
    }
    
    #convert into an average detection probability 
    #because no factors were found to affect detection probability, the value is constant
    for (j in 1:n.Transect){
    ESW.j[j] <- max(0.0001,ESW.jt[j,])
    }
    
    ESW.constant <- max(ESW.j)
    averagePa <- ESW.constant/W

    }
    ",fill = TRUE)
sink()

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("intercept.lt","averagePa","beta.forest","beta.military","beta.villages","beta.fields",
            "Density","nHat")

inits <- function(){list(sigma=runif(1,80,150))}

ni<-10000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
traceplot(out1)

library(ggmcmc)   
ggs_traceplot(ggs(as.mcmc(out1)))

########################################################################################

#Plotting predictions just on line transects

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits[sites.lt]<-out1$mean$Density[sites.lt]
myGridDF3km$fits[!1:nrow(myGridDF3km)%in%sites.lt]<-NA
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_lineonly.png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#get total number of predicted deer
sum(getValues(predRaster),na.rm=T)#1490.399

#########################################################################################


#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline.RData")

#########################################################################################