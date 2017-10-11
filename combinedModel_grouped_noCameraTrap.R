#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

######################################################################################

#######################
#Compile data for model#
########################

bugs.data<-list(#camera trap data
                #line transect data
                site.lt = sites.lt,
                n.Transect = length(unique(datafile$Transect)), 
                n.Yrs = length(unique(datafile$T)),
                W = 250,
                n = totalsInfo,
                n.Detections = nrow(detectionInfo),
                n.detectionSites = length(unique(detectionInfo$Site)),
                d.Forest = detectionInfo$forestcover,
                d.Military = detectionInfo$military,
                y = detectionInfo$Distance,
                d.Groupsize = detectionInfo$GroupSize,
                d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectAreas = transectAreas,
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                #covariates
                Forest=forestcover,
                Forest2=forestcover2,
                Military=military,
                Villages=villages,
                Water=waterMatrix,
                Month=monthMatrix)

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
sink("combinedModel_grouped.txt")
cat("
    model {
    
    #Common state model on abundance across all sites

    #Priors
    beta.forest ~ dnorm(0,0.001)
    beta.forest2 ~ dnorm(0,0.001)
    beta.villages ~ dnorm(0,0.001)
    beta.military ~ dnorm(0,0.001)

    #Model of factors affecting abundance
    for (i in 1:n.sites) { #across all sites   
      abund.effect[i] <-  beta.forest * Forest[i] + beta.military * Military[i] + beta.forest2 * Forest2[i] + 
                          beta.villages * Villages[i]
    }

    #Line transect data:

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
    ESW.constant <- max(ESW.j[])
    averagePa <- ESW.constant/W

    # MODEL NUMBER OF INDIVIDUALS DETECTED

    # PRIORS
    intercept.lt ~ dnorm(0,0.001)
    
    #Random year effect (we are using line transect data from 4 years)
    sd.n.time ~ dunif(0,3)
    tau.n.time <- pow(sd.n.time,-2)
    for(t in 1:n.Yrs){
      random.n.time[t] ~ dnorm(0,tau.n.time)
    }

    #fit model
    for (j in 1:n.Transect){
      for (t in 1:n.Yrs){
        n[j,t] ~ dpois(nHat[j,t]*transectAreas[j]*averagePa) 
        log(nHat[j,t]) <- intercept.lt + abund.effect[site.lt[j]] + random.n.time[t]
    }
    }

    #predict density (number of indivdiuals per grid cell) for all sites 
    #using the abundance effect and average group size
    #and the intercept of the line transect model for number of groups
    for(i in 1:n.sites){
      Density[i] <- exp(abund.effect[i] + intercept.lt)
    }

    #get predicted average density across whole area
    avDensity <- mean(Density)

    }
    ",fill = TRUE)
sink()

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

params <- c("b.gs.0","averagePa","beta.forest","beta.military","beta.forest2","beta.villages","avDensity","Density")

inits <- function(){list(z.ct = zst.ct,
                         a = ast,
                         sigma=runif(1,80,150))}

n.iter<-10000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
traceplot(out1)

########################################################################################

#Plotting predictions

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits<-out1$mean$Density
myGridDF$fits<-myGridDF3km$fits[match(myGridDF$Grid3km,myGridDF3km$Grid3km)]
out<-subset(myGridDF,!is.na(fits))
summary(out$fits)
predRaster<-rasterFromXYZ(out[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_noCameraTrap.png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#get total number of predicted deer
out2<-subset(out,!duplicated(Grid3km))
sum(out2$fits)#1159.388

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_noCameraTrap.RData")

#########################################################################################

