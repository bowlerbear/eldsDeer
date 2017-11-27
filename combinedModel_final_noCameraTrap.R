#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

######################################################################################
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')

#Use the JAGAM function to get the BUGS code for the spline
#http://www.petrkeil.com/?p=2385
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
library(mgcv)
jags.ready <- jagam(Deer.lt~1+s(x, y), 
                    data=subset(myGridDF3km,!is.na(Deer.lt)), 
                    family="poisson", 
                    sp.prior="log.uniform",
                    file="jagam.lt.txt")

#get the data bits we need from jags data
X = jags.ready$jags.data$X[,-1]
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero[-1]


#fit model
myGridDF3km$Deer.lt[sites.lt]
myGridDF3km$DeerCounts.lt<-NA
myGridDF3km$DeerCounts.lt[sites.lt]<-round(apply(groupInfo,1,mean))

model1 <- gam(DeerCounts.lt~1+s(x, y), 
                    data=subset(myGridDF3km,!is.na(Deer.lt)), 
                    family="poisson")

#smoothing terms are significant...

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
                y = detectionInfo$Distance,
                d.Groupsize = detectionInfo$GroupSize,
                #d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectAreas = transectAreas,
                #covariates
                Forest=as.numeric(scale(log(myGridDF3km$ForestCover+1)))[sites.lt],
                ForestB = forestcoverB,
                #Forest2=forestcover2,
                Fields=as.numeric(scale(sqrt(myGridDF3km$Fields+1)))[sites.lt],
                Military=myGridDF3km$MilitaryF[sites.lt],
                Villages=as.numeric(scale(sqrt(myGridDF3km$Villages+1)))[sites.lt],
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
    ## prior for s(x,y)... 
    K1 <- S1[1:29,1:29] * lambda[1]  + S1[1:29,30:58] * lambda[2]
    b[1:29] ~ dmnorm(zero[1:29],K1) 
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
      rho[i] ~ dunif(-1,1)
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
        nHat[j,t] <- Density[j,t]*transectAreas[j]*averagePa

        #(1)with covariates
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * Military[j] +
        #                  beta.fields * Fields[j] + beta.villages * Villages[j] 

        #(2)just random effects
        log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + 
                          random.s.year[t] + random.s.site[j]

        #(3)with spline
        #log(Density[j,t]) <- intercept.lt + eta[j]
        
        #(4) with spline and covariates
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * Military[j] +
        #                      beta.fields * Fields[j] + beta.villages * Villages[j] + eta[j]
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

    pi <- 3.141593
    
    #Priors
    sigma ~ dunif(50,140)
    b.d.0 ~ dunif(0,20)

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

params <- c("intercept.lt","intercept.groupsize",
            "averagePa","beta.forest","beta.military",
            "beta.fields","beta.villages",
            "site.s.sd","year.s.sd",
            "totDensity","n.transect","nHat",
            "fit","fit.new")

inits <- function(){list(sigma=runif(1,80,150))}

ni<-50000
nb<-10000
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
png(file="combinedModel_grouped_spline_noCameraTrap(2).png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline_noCameraTrap(2).RData")

#compare n and nHat
nMelt<-melt(groupInfo)
nHatMelt<-melt(out1$mean$nHat)
qplot(nMelt$value,nHatMelt$value)

#save similations
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#1 - 0.0023383
#2 - 0.2778167
#3 - 0.3454833
#4 - 0.0021833

#########################################################################################