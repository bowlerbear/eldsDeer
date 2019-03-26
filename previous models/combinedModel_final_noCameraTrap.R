#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

######################################################################################

##############################################
#Sort out the code for the spline model#
##############################################

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
#Use the JAGAM function to get the BUGS code for the spline

library(mgcv)
summary(m.gam<-gam(Deer~ 1 + s(x, y),data=myGridDF3km, family=binomial(link="cloglog")))
gam.check(m.gam)

jags.ready <- jagam(Deer~ 1 + s(x, y,bs="tp",k=10), 
                    data=myGridDF3km, 
                    family=binomial, 
                    sp.prior="log.uniform",
                    file="jagam.txt")


#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero


#######################
#Compile data for model#
########################

bugs.data<-list(#line transect data
  site.lt = sites.lt,
  n.Transect = length(unique(datafile$Transect)), 
  n.Yrs = length(unique(datafile$T)),
  W = 250,
  n = groupInfo,
  n.Detections = nrow(detectionInfo),
  #n.detectionSites = length(unique(detectionInfo$Site)),
  #d.Forest = as.numeric(scale(log(detectionInfo$forestcover+1))),
  #d.Military = ifelse(detectionInfo$military>0,1,0),
  d.transectId = as.numeric(factor(
    interaction(detectionInfo$transectID,detectionInfo$T))),
  n.d.transectId = length(unique(factor(
    interaction(detectionInfo$transectID,detectionInfo$T)))),
  y = detectionInfo$Distance,
  d.Groupsize = detectionInfo$GroupSize,
  #d.Site = detectionInfo$Site,
  n.TransectYrs = nrow(transectInfo),
  ty.combos = transectInfo,
  TransYrIdx = TransYrIdx,
  zeros.dist = rep(0,nrow(detectionInfo)),
  transectLengths = mytransectLengths,
  n.transectId = length(unique(myGridDF3km$transectId[!is.na(myGridDF3km$Deer.lt)])),
  transectId = myGridDF3km$transectId[!is.na(myGridDF3km$Deer.lt)],
  #total number of sites  
  n.sites = length(unique(myGridDF3km$Grid3km)),
  #covariates
  Forest=as.numeric(scale(log(myGridDF3km$ForestCover+1))),
  Forest2=as.numeric(scale(myGridDF3km$ForestCover^2)),
  ForestF=myGridDF3km$ForestCoverF,
  ForestF2=myGridDF3km$ForestCoverF2,
  Fields=as.numeric(scale(log(myGridDF3km$Fields+1))),
  Military=myGridDF3km$MilitaryF,
  MilitaryN=as.numeric(scale(myGridDF3km$Military)),
  Villages=as.numeric(scale(log(myGridDF3km$Villages+1))),
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
    
    #effect of covariates
    #gs.forest ~ dnorm(0,0.001)#no effect
    #gs.military ~ dnorm(0,0.001)#no effect
    
    
    for(i in 1:n.Detections){
    d.Groupsize[i] ~ dpois(exp.GroupSize[i])
    #log(exp.GroupSize[i]) <- intercept.groupsize + gs.forest * d.Forest[i] + 
    #                          gs.military * d.Military[i] 
    
    log(exp.GroupSize[i]) <- intercept.groupsize 
    }

    gp<-exp(intercept.groupsize)

    #Model of factors affecting abundance

    beta.forest2 ~ dnorm(0,0.001) 
    beta.forest ~ dnorm(0,0.001)  
    beta.military ~ dnorm(0,0.001)
    beta.fields ~ dnorm(0,0.001)
    beta.villages ~ dnorm(0,0.001)
    
    # Line transect data:

    # MODEL NUMBER OF GROUPS DETECTED

    # PRIORS
    intercept.lt ~ dnorm(0,0.001)

    #if including a spline?

    #Model for the spline
      #the linear predictor
    eta <- X %*% b ## linear predictor for the spline
    
    ## Parametric effect priors
    for (i in 1:1) { b[i] ~ dnorm(0,0.27) }
    
    ## prior for s(x,y)... 
    K1 <- S1[1:9,1:9] * lambda[1]  + S1[1:9,10:18] * lambda[2]
    b[2:10] ~ dmnorm(zero[2:10],K1) 
    
    ## smoothing parameter priors 
    for (i in 1:2) {
    rho[i] ~ dunif(-3,3)
    lambda[i] <- exp(rho[i])
    }


    #fit model
    for (j in 1:n.Transect){
      for (t in 1:n.Yrs){

        n[j,t] ~ dpois(nHat[j,t]) 
        nHat[j,t] <- (Density[j,t])*surveyArea[j,t]
        surveyArea[j,t]<-(transectLengths[j,t]*(ESW.constant/1000)*2)/9
        
        #(1)with all covariates
        #log(Density[j,t]) <- intercept.lt + beta.forest * Forest[j] + beta.military * MilitaryN[j] +
        #                  beta.fields * Fields[j] + beta.villages * Villages[j] 

        #(3)with spline
        #log(Density[j,t]) <- eta[j]
        
        #(4) with spline and covariates
        log(Density[j,t]) <- beta.forest * Forest[j] + beta.military * MilitaryN[j] +eta[j]
        
      }
    }

    #derived parameters
    for(j in 1:n.Transect){
      #get predicted total across whole area
      n.transect[j] <-mean(Density[j,]*gp)#mean across years 
    }
      totDensity <- sum(n.transect)


    #calculate the Bayesian p-value
    e <- 0.00000001
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
    b.group.size ~ dnorm(0,0.001)

   #random transect effect
    for(i in 1:n.d.transectId){
        random.transect[i] ~ dnorm(0,random.transect.tau)
    }
    random.transect.tau <- pow(random.transect.sd,-2)
    random.transect.sd ~ dunif(0,10)

    ##### Begin model for *all detections*
    
    for( i in 1:n.Detections){
    
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    mu.df[i] <- b.d.0
    #mu.df[i] <- b.d.0 + b.group.size * d.Groupsize[i]
    #mu.df[i] <- b.d.0 + random.transect[d.transectId[i]]

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

    #predict over total area
    for(j in 1:n.sites){
      log(predDensity[j]) <- intercept.lt + beta.forest * Forest_All[j] + 
                              beta.military * Military_All[j]
      
      totalIndiv[j] <- predDensity[j]*gp
    }

    totalPredDensity <- sum(totalIndiv)
    }

    ",fill = TRUE)
sink()

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

params <- c("ESW.constant",
            "rho","intercept.lt","intercept.groupsize",
            "averagePa","beta.forest","beta.military",
            "beta.fields","beta.villages",
            "site.s.sd","year.s.sd",
            "totalPredDensity",
            "totDensity","n.transect","nHat",
            "fit","fit.new","totalIndiv")

inits <- function(){list(sigma=runif(1,80,150))}

ni<-50000
nb<-10000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
save(out1,file="out1_final_noCameraTrap(5).RData")

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
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,170))
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
#3 - 0.4144617
#4 - 0.00218
#5 - 0.05243111 
#6 - 0.002555556
#7 - 0.03

#########################################################################################

#for predictions across the whole area
myGridDF3km$fits<-out1$mean$totalIndiv
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_noCameraTrap(5).png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,160))
plot(sws,add=T)
dev.off()

#########################################################################################