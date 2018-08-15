#################################################################################

###################
#Retrive the data##
###################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

####################################################################################################################

##############################################
#Sort out the code for the spline model#
##############################################

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
#Use the JAGAM function to get the BUGS code for the spline

library(mgcv)
summary(m.gam<-gam(Deer~ 1 + s(x, y),data=myGridDF3km, family=binomial(link="cloglog")))
gam.check(m.gam)

jags.ready <- jagam(Deer~ 1 + s(x, y), 
                    data=myGridDF3km, 
                    family=binomial, 
                    sp.prior="log.uniform",
                    file="jagam.txt")


#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero

##################################################################################

#add NAs to the groupInfo data

new.my.n <-matrix(data=NA,nrow=nrow(myGridDF3km),ncol=3)
new.my.n[sites.lt,]<-groupInfo
new.transectlengths <-matrix(data=NA,nrow=nrow(myGridDF3km),ncol=3)
new.transectlengths[sites.lt,]<-mytransectLengths
new.transectlengths[is.na(new.transectlengths)]<-0


#index all data points
transectInfo<-data.frame(expand.grid(Transect=1:nrow(myGridDF3km),T=1:3))
TransYrIdx<-matrix(nrow=nrow(datafileObs),ncol=nrow(transectInfo))
TransYrIdx[]<-0
for(i in 1:nrow(datafileObs)){
  TransYrIdx[i,which(transectInfo$T==datafileObs$T[i]&
                       datafile$Transect[i]==transectInfo$Transect)]<-1
}

#####################################
#Compile the remaining data for model#
#####################################

bugs.data<-list(#line transect data
                site.lt = sites.lt,
                n.Transect = length(unique(datafile$Transect)), 
                n.Yrs = length(unique(datafile$T)),
                W = 250,
                n = new.my.n,
                n.Detections = nrow(detectionInfo),
                #n.detectionSites = length(unique(detectionInfo$Site)),
                d.Forest = as.numeric(scale(detectionInfo$forestcover)),
                d.Military = as.numeric(scale(detectionInfo$military)),
                d.transectId = as.numeric(factor(
                  interaction(detectionInfo$transectID,detectionInfo$T))),
                n.d.transectId = length(unique(factor(
                  interaction(detectionInfo$transectID,detectionInfo$T)))),
                y = detectionInfo$Distance,
                d.Groupsize = detectionInfo$GroupSize,
                d.GroupsizeS = as.numeric(detectionInfo$GroupSize-median(detectionInfo$GroupSize)),
                #d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectLengths = new.transectlengths,
                n.transectId = length(unique(myGridDF3km$transectId[!is.na(myGridDF3km$Deer.lt)])),
                transectId = myGridDF3km$transectId[!is.na(myGridDF3km$Deer.lt)],
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                #covariates
                Forest=as.numeric(scale((myGridDF3km$ForestCover))),
                Forest2=as.numeric(scale(log(myGridDF3km$ForestCover+1)^2)),
                ForestF=myGridDF3km$ForestCoverF,
                ForestF2=myGridDF3km$ForestCoverF2,
                Fields=as.numeric(scale(log(myGridDF3km$Fields+1))),
                Military=myGridDF3km$MilitaryF,
                MilitaryN=as.numeric(scale(myGridDF3km$Military)),
                Villages=as.numeric(scale(log(myGridDF3km$Villages+1))),
                Water=waterMatrix,
                Lure=lureMatrix,
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
    beta.forest2 ~ dnorm(0,0.001)     
    beta.military ~ dnorm(0,0.001)
    beta.fields ~ dnorm(0,0.001)
    beta.villages ~ dnorm(0,0.001)

     #Model for the spline
      #the linear predictor
      eta <- X %*% b ## linear predictor for the spline
  
      ## Parametric effect priors
      for (i in 1:1) { b[i] ~ dnorm(0,0.01) }
       
      ## prior for s(x,y)... 
      K1 <- S1[1:29,1:29] * lambda[1]  + S1[1:29,30:58] * lambda[2]
      b[2:30] ~ dmnorm(zero[2:30],K1) 

      ## smoothing parameter priors 
      for (i in 1:2) {
        rho[i] ~ dunif(-1,1)
        lambda[i] <- exp(rho[i])
      }
  
    
    #(2) Line transect data:

    # MODEL GROUP SIZE

    intercept.groupsize ~ dnorm(0,0.001)
    
    for(i in 1:n.Detections){
      d.Groupsize[i] ~ dpois(exp.GroupSize[i])
      log(exp.GroupSize[i]) <- intercept.groupsize 
    }

    # MODEL NUMBER OF GROUPS DETECTED

    # PRIORS
    intercept.lt ~ dunif(0,10)

    #fit model
    for (j in 1:n.sites){
      for (t in 1:n.Yrs){

        n[j,t] ~ dpois(nHat[j,t]) 

        #work out area of each cell that is surveyed
        surveyArea[j,t]<-((transectLengths[j,t]/1000)*(ESW[j,t]/1000)*2)/9

        #relate fraction seen to fraction available in whole cell
        nHat[j,t] <- (expN[j,t])*surveyArea[j,t]

        log(expN[j,t]) <- intercept.lt + beta.forest * Forest[j] +
                          beta.military * MilitaryN[j] +
                          eta[j] + beta.fields* Fields[j]

        #log(expN[j,t]) <- intercept.lt + eta[j]

      }
    }

    #predict density (number of indivdiuals per grid cell) for all sites 
    #using the abundance effect and average group size
    #and the intercept of the line transect model for number of groups
    for(i in 1:n.sites){
      Density[i] <- mean(expN[i,]) * exp(intercept.groupsize) 
    }

    #get predicted total across whole area
    avDensity <- mean(Density)
    totDensity <- sum(Density)

    #calculate the Bayesian p-value
    e <- 0.0001
    for(j in 1:n.sites){
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
    sigma ~ dunif(10,200)
    b.d.0 ~ dunif(0,6)

    #random transect effect
    for(i in 1:n.d.transectId){
        random.transect[i] ~ dnorm(0,random.transect.tau)
    }
    random.transect.tau <- pow(random.transect.sd,-2)
    random.transect.sd ~ dunif(0,5)

    ##### Begin model for *all detections*
    
    for(i in 1:n.Detections){
    
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    mu.df[i] <- b.d.0 + random.transect[d.transectId[i]] 

    # estimate of sd and var, given coeffs above
    sig.df[i] <- exp(mu.df[i])
    sig2.df[i] <- sig.df[i]*sig.df[i]
    
    # effective strip width
    esw[i] <- sqrt(pi * sig2.df[i] / 2)
    f0[i] <- 1/esw[i] #assuming all detected on the line
    
    # LIKELIHOOD
    # estimate sigma using zeros trick
    #y[i] ~ dunif(0,W)
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
    
    #convert into detection probability 
    mean.esw <- sqrt(pi * pow(exp(b.d.0),2) / 2) 
    for(j in 1:n.sites){
      for(t in 1:n.Yrs){
        ESW[j,t] <- ifelse(equals(ESW.jt[j,t],0),mean.esw,ESW.jt[j,t])
      }
    }

    }
    ",fill = TRUE)
sink()

###############
#Run the model#
###############

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("random.transect","b.d.0",
            "b.group.size","mean.esw",
            "random.transect.sd","intercept.year",
            "intercept.lt","intercept.ct","intercept.groupsize",
            "averagePa","beta.forest","beta.forest2","beta.military",
            "abund.slope",
            "beta.villages","beta.fields",
            "alpha.month","theta.water","theta.lure",
            "gs.forest","gs.military",
            "avDensity","totDensity","Density",
            "nHat","fit","fit.new")

zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

inits <- function(){list(sigma=runif(1,80,150))}

ni<-10000
nb<-3000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni,parallel=TRUE)

print(out1,2)
traceplot(out1)
save(out1,file="out1_final_ltonly.RData")

########################################################################################

#Plotting predictions
library(gplot)
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits<-out1$mean$Density
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_ltonly.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,115))
plot(sws,add=T)
dev.off()

#spline onle
library(gplot)
setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits<-NA
myGridDF3km$fits[sites.lt]<-out1$mean$Density[sites.lt]
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_splineonly_ltonly.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,115))
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline.RData")

#save similations for Bayesian p-value
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#0.3349524, ForestF
#0.4086667, log Forest
#0.3629524, Forest
#########################################################################################