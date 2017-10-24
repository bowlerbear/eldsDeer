#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

#####################################
#Compile the remaining data for model#
#####################################

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
                Forest=forestcoverB,
                #Forest2=forestcover2,
                #Fields=as.numeric(scale(sqrt(fields+1))),
                Military=ifelse(military>0,1,0),
                #Villages=as.numeric(scale(sqrt(villages+1))),
                Water=waterMatrix,
                Month=monthMatrix)

#################
#Write the model#
#################

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
sink("combinedModel_grouped_spline.txt")
cat("
    model {
    
    #Model of factors affecting abundance
    beta.forest ~ dnorm(0,0.001)  
    beta.military ~ dnorm(0,0.001)
    #beta.fields ~ dnorm(0,0.001)
    #beta.villages ~ dnorm(0,0.001)
    
    for (i in 1:n.sites) { #across all sites   
      abund.effect[i] <-  beta.forest * Forest[i] + beta.military * Military[i]
    }

    #Different observation models depending on the data:

    #(1) Camera trap data
    
    #Priors
    mean.p.ct ~ dunif(0, 1)
    alpha.ct <- cloglog(mean.p.ct)
    theta.water ~ dnorm(0,0.001)
    alpha.month ~ dnorm(0,0.001)
    intercept.ct ~ dnorm(0,0.001)

    # Intercepts availability probability
    for(j in 1:n.1kmGrids){
    int.theta[j] ~ dunif(0,1) 
    }
    
    #Model (multi-scale occupancy model using cloglog so we are esssentially modelling abundance!)
    for (i in 1:n.CameraTrapSites) { 
    z.ct[i] ~ dbern(psi[i])
    cloglog(psi[i]) <- intercept.ct + abund.effect[site.ct[i]]
    
    #availbility for detection model
    for (j in 1:n.1kmGrids){
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z.ct[i] * theta[i,j]
    cloglog(theta[i,j]) <- logit(int.theta[j]) + theta.water * Water[i,j]
    
    #detection submodel 
    for (k in 1:n.reps) { # Loop over replicate surveys
    y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
    mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
    cloglog(p.ct[i,j,k]) <- alpha.ct + alpha.month * Month[i,j,k] + abund.effect[site.ct[i]]
    }
    }
    }

    #predict density (number of indivdiuals per grid cell) for all sites 
    #using the abundance effect and average group size
    #and the intercept of the line transect model for number of groups
    for(i in 1:n.sites){
    cloglog(Density[i]) <- intercept.ct + abund.effect[i]
    }
    
    #get predicted total across whole area
    avDensity <- mean(Density)
    totDensity <- sum(Density)
    

    }
    ",fill = TRUE)
sink()

###############
#Run the model#
###############

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("intercept.ct","mean.p.ct",
            "beta.forest","beta.military",
            #"beta.villages","beta.fields",
            "avDensity","totDensity","Density")

inits <- function(){list(intercept.ct=runif(1,8,15))}

ni<-20000
nb<-5000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
traceplot(out1)

library(ggmcmc)   
ggs_traceplot(ggs(as.mcmc(out1)))

########################################################################################

#Plotting predictions

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits<-out1$mean$Density
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_noLineTransect.png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline_noLineTransect.RData")

#########################################################################################