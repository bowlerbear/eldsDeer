#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

##############################################
#First sort out the code for the spline model#
##############################################

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')

#Use the JAGAM function to get the BUGS code for the spline
#http://www.petrkeil.com/?p=2385
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
library(mgcv)
#myGridDF3km$ForestS<-as.numeric(scale(log(forestcover+1)))

jags.ready <- jagam(Deer~ 1 + s(x, y), 
                    data=myGridDF3km, 
                    family="poisson", 
                    sp.prior="log.uniform",
                    file="jagam.txt")


#get the data bits we need from jags data
X = jags.ready$jags.data$X[,-1]
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero[-1]

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
                n = groupInfo,
                n.Detections = nrow(detectionInfo),
                #n.detectionSites = length(unique(detectionInfo$Site)),
                #d.Forest = as.numeric(scale(log(detectionInfo$forestcover+1))),
                #d.Military = ifelse(detectionInfo$military>0,1,0),
                y = detectionInfo$Distance,
                d.Groupsize = detectionInfo$GroupSize,
                #d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectAreas = transectAreas,
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                #covariates
                Forest=as.numeric(scale(log(myGridDF3km$ForestCover+1))),
                #Forest2=forestcover2,
                Fields=as.numeric(scale(sqrt(myGridDF3km$Fields+1))),
                Military=myGridDF3km$MilitaryF,
                Villages=as.numeric(scale(sqrt(myGridDF3km$Villages+1))),
                Water=waterMatrix,
                Month=monthMatrix,
                #spatial covariates
                X = X,
                S1 = S1,
                zero = zero)

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
    beta.fields ~ dnorm(0,0.001)
    beta.villages ~ dnorm(0,0.001)

    #Common model
    for (i in 1:n.sites) { #across all sites   
      #(1)
      #abund.effect[i] <-  eta[i] + beta.forest * Forest[i] + beta.military * Military[i]
      #(2)
      #abund.effect[i] <-  beta.forest * Forest[i] + beta.military * Military[i]
      #(3)
      #abund.effect[i] <-  beta.forest * Forest[i] + beta.military * Military[i] + 
      #                    beta.fields * Fields[i] + beta.villages * Villages[i]
      #(4)
      #abund.effect[i] <-  beta.forest * Forest[i] + beta.military * Military[i] + 
      #                    beta.fields * Fields[i] + beta.villages * Villages[i] +
      #                    eta[i]
      #(5)
      #abund.effect[i] <-  beta.forest * Forest[i] + beta.military * Military[i] + 
      #                     beta.villages * Villages[i] + eta[i]
      #(6)
      abund.effect[i] <-  eta[i]
      }

    #Model for the spline
    #the linear predictor
    eta <- X %*% b ## linear predictor for the spline

    ## prior for s(x,y)... 
    K1 <- S1[1:29,1:29] * lambda[1]  + S1[1:29,30:58] * lambda[2]
    b[1:29] ~ dmnorm(zero[1:29],K1) 
    
    ## smoothing parameter priors.
    for (i in 1:2) {
      rho[i] ~ dunif(-1,1)
      lambda[i] <- exp(rho[i])
    }

    #Different observation models depending on the data:

    #(1) Camera trap data
    
    #Priors
    theta.water ~ dnorm(0,0.001)
    alpha.month ~ dnorm(0,0.001)
    abund.slope ~ dnorm(0,0.001)
    abund.slope2 ~ dnorm(0,0.001)

    #intercepts for each model
    mean.p.ct ~ dunif(0,1)
    intercept.ct ~ dunif(0,1)
    int.theta ~ dunif(0,1) 


    #Model (multi-scale occupancy model using cloglog so we are esssentially modelling abundance!)
    for (i in 1:n.CameraTrapSites) { 
    z.ct[i] ~ dbern(psi[i])
    cloglog(psi[i]) <- cloglog(intercept.ct) + abund.effect[site.ct[i]]
    
    #availbility for detection model
    for (j in 1:n.1kmGrids){
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z.ct[i] * theta[i,j]
    #cloglog(theta[i,j]) <- cloglog(int.theta) + theta.water * Water[i,j]
    cloglog(theta[i,j]) <- cloglog(int.theta)

    #detection submodel 
    for (k in 1:n.reps) { # Loop over replicate surveys
    y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
    mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
    #cloglog(p.ct[i,j,k]) <- cloglog(mean.p.ct) + abund.slope * log(Density[site.ct[i]]) + alpha.month * Month[i,j,k]
    cloglog(p.ct[i,j,k]) <- cloglog(mean.p.ct) + abund.slope * log(Density[site.ct[i]])
    
    }
    }
    }
    
    #(2) Line transect data:

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

    # MODEL NUMBER OF GROUPS DETECTED

    # PRIORS
    intercept.lt ~ dunif(0,10)

    #fit model
    for (j in 1:n.Transect){
      for (t in 1:n.Yrs){
        n[j,t] ~ dpois(nHat[j,t]) 
        nHat[j,t] <- expN[j,t]*transectAreas[j]*averagePa
        log(expN[j,t]) <- intercept.lt + abund.effect[site.lt[j]]
    }
    }

    #predict density (number of indivdiuals per grid cell) for all sites 
    #using the abundance effect and average group size
    #and the intercept of the line transect model for number of groups
    for(i in 1:n.sites){
      log(Density[i]) <- (intercept.lt + abund.effect[i])*(intercept.groupsize) 
    }

    #get predicted total across whole area
    avDensity <- mean(Density)
    totDensity <- sum(Density)



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

###############
#Run the model#
###############

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("intercept.lt","intercept.ct","lambda","intercept.groupsize",
            "rho","averagePa","beta.forest","beta.military","abund.slope","abund.slope2",
            "beta.villages","beta.fields",
            "alpha.month","theta.water","obs.sd",
            "gs.forest","gs.military",
            "avDensity","totDensity","Density","nHat","fit","fit.new")

zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

inits <- function(){list(sigma=runif(1,80,150),
                         z.ct = zst.ct, 
                         a = ast)}

ni<-50000
nb<-10000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
traceplot(out1)

library(ggmcmc)   
out2<-ggs(out1$samples)
ggs_traceplot(filter(out2,Parameter%in%c("intercept.lt","intercept.ct","totDensity")))

########################################################################################

#Plotting predictions

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits<-out1$mean$Density
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline(6).png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline(6).RData")

#save similations
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)

#3, 0.00205
#4, 0.163883
#6, 0.3200167

#########################################################################################