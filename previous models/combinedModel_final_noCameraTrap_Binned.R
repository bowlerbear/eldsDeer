#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

######################################################################################

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')

#Use the JAGAM function to get the BUGS code for the spline
#http://www.petrkeil.com/?p=2385
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
library(mgcv)
jags.ready <- jagam(Deer.lt~1+s(x, y,k=10), 
                    data=subset(myGridDF3km,!is.na(Deer.lt)), 
                    family="poisson", 
                    sp.prior="log.uniform",
                    file="jagam.lt.txt")

#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero

#fit model
myGridDF3km$Deer.lt[sites.lt]
myGridDF3km$DeerCounts.lt<-NA
myGridDF3km$DeerCounts.lt[sites.lt]<-round(apply(groupInfo,1,mean))
model1 <- gam(DeerCounts.lt~1+s(x, y,k=10), 
                    data=subset(myGridDF3km,!is.na(Deer.lt)), 
                    family="poisson")
summary(model1)
#smoothing terms are significant...

######################################################################################

#bin distance data
summary(detectionInfo$Distance)
hist(detectionInfo$Distance)
Breaks=c(0,50,100,150,200,250)
Midpt= (Breaks[2:6]+Breaks[1:5])/2
detectionInfo$Distance<-cut(detectionInfo$Distance,breaks=Breaks)
levels(detectionInfo$Distance) <-1:length(unique(detectionInfo$Distance))
detectionInfo$Distance<-as.numeric(as.character(detectionInfo$Distance))
detectionInfo$Distance[is.na(detectionInfo$Distance)]<-1
unique(detectionInfo$Distance)
hist(detectionInfo$Distance)

######################################################################################

#######################
#Compile data for model#
########################

bugs.data<-list(
                #line transect data
                site.lt = sites.lt,
                n.Transect = length(unique(datafile$Transect)), 
                n.Yrs = length(unique(datafile$T)),
                W = 250,
                n = groupInfo,
                n.Detections = nrow(detectionInfo),
                #n.detectionSites = length(unique(detectionInfo$Site)),
                #d.Forest = detectionInfo$forestcover,
                #d.Military = detectionInfo$military,
                dclass = detectionInfo$Distance,
                B=max(Breaks),
                nD =length(Breaks)-1,
                midpt=Midpt,
                delta=Breaks[2]-Breaks[1],
                d.Groupsize = detectionInfo$GroupSize,
                #d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectLengths = myGridDF3km[sites.lt,c("t2014","t2015","t2016")],
                #covariates
                Forest=as.numeric(scale(log(myGridDF3km$ForestCover+1)))[sites.lt],
                ForestB=myGridDF3km$ForestCoverF[sites.lt],
                #Forest2=forestcover2,
                Fields=as.numeric(scale(sqrt(myGridDF3km$Fields+1)))[sites.lt],
                Military=myGridDF3km$MilitaryF[sites.lt],
                Villages=as.numeric(scale(sqrt(myGridDF3km$Villages+1)))[sites.lt],
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                ForestB_All=myGridDF3km$ForestCoverF,
                Forest_All=as.numeric(scale(log(myGridDF3km$ForestCover+1))),
                Villages_All=as.numeric(scale(sqrt(myGridDF3km$Villages+1))),
                Military_All=myGridDF3km$MilitaryF,
                #spatial covariates
                X = X,
                S1 = S1,
                zero = zero)


setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
sink("combinedModel_grouped_spline.txt")
cat("
    model {
    
    # MODEL GROUP SIZE

    intercept.groupsize ~ dnorm(0,0.001)
    obs.sd ~ dunif(0,10)
    obs.tau <- pow(obs.sd,-2)
    for(i in 1:n.Detections){
    random.obs[i] ~ dnorm(0,obs.tau)
    }
    
    #effect of covariates
    #gs.forest ~ dnorm(0,0.001)#no effect
    #gs.military ~ dnorm(0,0.001)#no effect
    
    for(i in 1:n.Detections){
    d.Groupsize[i] ~ dpois(exp.GroupSize[i])
    #log(exp.GroupSize[i]) <- intercept.groupsize + gs.forest * d.Forest[i] + 
    #                          gs.military * d.Military[i] + random.obs[i] 
    
    log(exp.GroupSize[i]) <- intercept.groupsize
    }


    #Model of factors affecting abundance
    beta.forest ~ dnorm(0,0.001)  
    beta.military ~ dnorm(0,0.001)
    beta.fields ~ dnorm(0,0.001)
    beta.villages ~ dnorm(0,0.001)
    
    # Line transect data:

    # MODEL NUMBER OF GROUPS DETECTED

    # PRIORS
    intercept.lt ~ dnorm(0,0.001)

    #if including a spline?
    eta <- X %*% b 

    ## Parametric effect priors CHECK tau=1/10^2 is appropriate!
    for (i in 1:1) { b[i] ~ dnorm(0,0.01) }
    
    ## prior for s(x,y)... 
    K1 <- S1[1:9,1:9] * lambda[1]  + S1[1:9,10:18] * lambda[2]
    b[2:10] ~ dmnorm(zero[2:10],K1) 
    
    ## smoothing parameter priors
    for (i in 1:2) {
      rho[i] ~ dunif(-3,3)
      lambda[i] <- exp(rho[i])
    }

    #random year effect
    year.s.sd ~ dunif(0,10)
    year.s.tau <- pow(year.s.sd,-2)
    for(t in 1:n.Yrs){
      random.s.year[t] ~ dunif(0,year.s.tau)
    }

    #random site effect
    site.s.sd ~ dunif(0,10)
    site.s.tau <- pow(site.s.sd,-2)
    for(j in 1:n.Transect){
      random.s.site[j] ~ dunif(0,site.s.tau)
    }

    #fit model
    for (j in 1:n.Transect){
      for (t in 1:n.Yrs){
        n[j,t] ~ dpois(nHat[j,t]) 
        nHat[j,t] <- (Density[j,t])*surveyArea[j,t]
        surveyArea[j,t]<-(transectLengths[j,t]*(ESW.constant/1000)*2)/9
        
        #(1)with all covariates
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * Military[j] +
        #                  beta.fields * Fields[j] + beta.villages * Villages[j] 

        #(2)just random effects
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + 
        #                  random.s.year[t] + random.s.site[j]

        #(3)with spline
        #log(Density[j,t]) <- eta[j]
        
        #(4) with spline and covariates
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * Military[j] +
        #                      beta.villages * Villages[j]+eta[j]
        
        #(5)#with only significant covariates
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * Military[j] +
        #                      beta.villages * Villages[j]
        
        #(6)with binary forest cover
        #log(Density[j,t]) <- intercept.lt + beta.forest * ForestB[j] + beta.military * Military[j]

        #(7)with forest cover and military
        log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * Military[j]

      }
    }

    #derived parameters
    gp<-exp(intercept.groupsize)
    for(j in 1:n.Transect){
      #get predicted total across whole area
      n.transect[j] <-mean(Density[j,]*gp)#mean across years 
    }
      totDensity <- sum(n.transect)


    #calculate the Bayesian p-value
    e <- 0.0001
    for(j in 1:n.Transect){
      for(t in 1:n.Yrs){
    # Fit assessments: Chi-squared test statistic and posterior predictive check
    chi2[j,t] <- pow((n[j,t] - nHat[j,t]),2) / (sqrt(nHat[j,t])+e)         # obs.
    nHat.new[j,t] ~ dpois(nHat[j,t])      # Replicate (new) data set
    chi2.new[j,t] <- pow((nHat.new[j,t] - nHat[j,t]),2) / (sqrt(nHat[j,t])+e) # exp
      }
    }

    # Add up discrepancy measures for entire data set
    for(t in 1:n.Yrs){
      fit.t[t] <- sum(chi2[,t])                     
      fit.new.t[t] <- sum(chi2.new[,t])             
    }

    fit <- sum(fit.t[])                     # Omnibus test statistic actual data
    fit.new <- sum(fit.new.t[])             # Omnibus test statistic replicate data

    
# DETECTABILITY COMPONENT#####

    #Priors
    sigma ~ dunif(50,200)

    ##### Begin model for *all detections*
  
    for( i in 1:n.Detections){
      dclass[i] ~ dcat(fc)
    }


    # Compute detection probability
    for(k in 1:nD){#for each distance class
      log(p[k]) <- -midpt[k]*midpt[k]/(2*sigma*sigma)
      pi[k] <- delta/B
      f[k] <- p[k]*pi[k]
      fc[k] <- f[k]/pcap
    }

    pcap <- sum(f)  # Overall detection probability
    
  
    averagePa <-pcap
    ESW.constant <- averagePa*W
    

    #predict over total area
    #for(j in 1:n.sites){
    #  log(predDensity[j]) <- intercept.lt + beta.forest * Forest_All[j] + 
    #                          beta.military * Military_All[j] +
    #                          beta.villages * Villages_All[j]
    #  totalIndiv[j] <- predDensity[j]*gp
    #}

    #totalPredDensity <- sum(totalIndiv)
    
    }
    ",fill = TRUE)
sink()

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

params <- c("ESW.constant","rho","intercept.lt","intercept.groupsize",
            "averagePa","beta.forest","beta.military",
            "beta.fields","beta.villages",
            "site.s.sd","year.s.sd",
            "totalPredDensity",
            "totDensity","n.transect","nHat",
            "fit","fit.new","totalIndiv")

inits <- function(){list(sigma=runif(1,80,150))}

ni<-50000
nb<-5000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)

#Plots
traceplot(out1)
library(ggmcmc)   
ggs_traceplot(ggs(as.mcmc(out1)))

########################################################################################

#Plotting predictions just on line transects

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
#add fits
myGridDF3km$fits<-NA
myGridDF3km$fits[sites.lt]<-out1$mean$n.transect

#plot them
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_noCameraTrap(3).png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline_noCameraTrap(3).RData")

#compare n and nHat
nMelt<-melt(groupInfo)
nHatMelt<-melt(out1$mean$nHat)
qplot(nMelt$value,nHatMelt$value)

#save similations
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#1 - 0.04743
#2 - 0.13625
#3 - 0.46431
#4 - 0.00218
#5 - 0.05221
#6 - 0.002555556
#7 - 0.03

#########################################################################################

#for predictions across the whole area
myGridDF3km$fits<-out1$mean$totalIndiv
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_noCameraTrap(5).png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#########################################################################################