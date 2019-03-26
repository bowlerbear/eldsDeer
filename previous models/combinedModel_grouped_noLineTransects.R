#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

######################################################################################

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
                Forest=forestcoverB,
                Forest2=forestcover2,
                Fields=fields,
                Military=military,
                Villages=villages,
                Water=waterMatrix,
                Month=monthMatrix)

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
sink("combinedModel.txt")
cat("
    model {
    
    #Common state model on abundance across all sites

    #Priors
    beta.forest ~ dnorm(0,0.001)
    beta.forest2 ~ dnorm(0,0.001)
    beta.villages ~ dnorm(0,0.001)
    beta.military ~ dnorm(0,0.001)
    beta.fields ~ dnorm(0,0.001)


    #Model of factors affecting abundance
    for (i in 1:n.sites) { #across all sites   
      abund.effect[i] <-  beta.forest * Forest[i]+ beta.military * Military[i] + beta.fields * Fields[i] + beta.villages * Villages[i]
    }

    #Different observation models depending on the data:

    #(1) Camera trap data

    #Priors

    mean.p.ct ~ dunif(0, 1)
    alpha.ct <- logit(mean.p.ct)
    intercept.ct ~ dnorm(0,0.001)
    theta.water ~ dnorm(0,0.001)
    alpha.month ~ dnorm(0,0.001)

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
    cloglog(theta[i,j]) <- logit(int.theta[j]) + theta.water * Water[i,j]+ abund.effect[site.ct[i]]
  
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
      log(Density[i]) <- abund.effect[i] + intercept.ct
    }
    
    #get predicted average density across whole area
    #avDensity <- mean(Density)
    #totDensity <- sum(Density)

    }
    ",fill = TRUE)
sink()

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

params <- c("intercept.ct","averagePa","beta.forest","beta.military","beta.villages","beta.fields",
            "abund.effect","Density")

inits <- function(){list(z.ct = zst.ct,
                         a = ast)}

ni<-10000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel.txt", n.thin=nt,
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
png(file="combinedModel_grouped_noLineTransects.png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#get total number of predicted deer
out2<-subset(out,!duplicated(Grid3km))
sum(out2$fits)

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_noLineTransects.RData")

#########################################################################################

#with all cloglog and abund effect on detection probability
                mean    sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
beta.forest    -3.01  0.67  -4.35  -3.00  -1.66    FALSE 1.00 1.09    29
beta.military   0.19  0.06   0.07   0.19   0.31    FALSE 1.00 1.04    58
beta.villages  -9.13  3.09 -15.33  -9.08  -3.34    FALSE 1.00 1.01   301
beta.fields    -0.13  0.28  -0.66  -0.14   0.45     TRUE 0.69 1.05    52

#with all cloglog and abund effect on the state model
                mean    sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
beta.forest   -13.83 12.07 -43.16  -7.10  -2.38    FALSE 1.00 2.98     4
beta.military  24.38 18.68   0.70  20.55  67.90    FALSE 0.99 1.01   358
beta.villages -11.74 19.50 -56.45  -9.24  21.53     TRUE 0.72 1.16    20
beta.fields    -6.75  5.85 -21.12  -4.74  -0.21    FALSE 0.99 2.24     5

#with all cloglog and abund effect on all
                mean   sd   2.5%    50%  97.5% overlap0    f Rhat n.eff
beta.forest    -1.81 0.45  -2.78  -1.78  -1.01    FALSE 1.00 1.00   910
beta.military   0.17 0.05   0.07   0.17   0.27    FALSE 1.00 1.00   656
beta.villages  -6.47 2.27 -10.99  -6.43  -2.12    FALSE 1.00 1.00   766
beta.fields    -0.27 0.19  -0.65  -0.26   0.07     TRUE 0.93 1.00  1059

#######################################################################################

