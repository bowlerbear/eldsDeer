#################################################################################

###################
#Retrive the data##
###################

source('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')
old <- myGridDF3km
####################################################################################################################

##############################################
#Sort out the code for the spline model#
##############################################
myGridDF3km <- old

setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')
#Use the JAGAM function to get the BUGS code for the spline
myGridDF3km$Deer.lt <- NA
myGridDF3km$Deer.lt[sites.lt] <- round(apply(groupInfo,1,mean))
medGridDF3km.x <- median(myGridDF3km$x[!is.na(myGridDF3km$Deer.lt)])
medGridDF3km.y <- median(myGridDF3km$y[!is.na(myGridDF3km$Deer.lt)])
myGridDF3km$x <- myGridDF3km$x - medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y - medGridDF3km.y

library(mgcv)
df <- subset(myGridDF3km,!is.na(Deer.lt))
m.gam<-gam(Deer.lt ~ s(x, y,k=4),data=df, family=poisson)
gam.check(m.gam)
summary(m.gam)
#spline to residuals
df$resids <- residuals(m.gam)
m.gam<-gam(resids ~ s(x, y),data=df)
summary(m.gam)

jags.ready <- jagam(Deer.lt~ 1 + s(x, y,k=4), 
                    data=subset(myGridDF3km,!is.na(Deer.lt)), 
                    family=poisson, 
                    sp.prior="log.uniform",
                    file="jagam.txt")


#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero

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
  n.detectionSites = length(unique(detectionInfo$Site)),
  detectionSites = detectionInfo$Transect,
  d.Forest = as.numeric(scale(sqrt(detectionInfo$forestcover+0.01))),
  d.Military = as.numeric(scale(sqrt(detectionInfo$militaryN+0.01))),
  d.Villages = as.numeric(scale(sqrt(detectionInfo$village+0.01))),
  d.Field = as.numeric(scale(sqrt(detectionInfo$field+0.01))),
  d.transectId = as.numeric(factor(
    interaction(detectionInfo$transectID,detectionInfo$T))),
  n.d.transectId = length(unique(factor(
    interaction(detectionInfo$transectID,detectionInfo$T)))),
  y = detectionInfo$Distance,
  d.Groupsize = detectionInfo$GroupSize,
  d.GroupsizeS = as.numeric(detectionInfo$GroupSize-median(detectionInfo$GroupSize)),
  d.Year = detectionInfo$T,
  d.Site = detectionInfo$Site,
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
  Forest=as.numeric(scale(sqrt(myGridDF3km$ForestCover+0.01))),
  Forest2=as.numeric(scale(sqrt(myGridDF3km$ForestCover^2+0.01))),
  ForestF=myGridDF3km$ForestCoverF,
  ForestF2=myGridDF3km$ForestCoverF2,
  Fields=as.numeric(scale(sqrt(myGridDF3km$Fields+0.01))),
  Military=myGridDF3km$MilitaryF,
  MilitaryN=as.numeric(scale(sqrt(myGridDF3km$Military+0.01))),
  Villages=as.numeric(scale(sqrt(myGridDF3km$Village+0.01))),
  Water=waterMatrix,
  Year.mat=yearMatrix,
  Lure=lureMatrix,
  Month=monthMatrix,
  #spatial covariates
  X = X,
  S1 = S1,
  zero = zero,
  nspline = length(zero),
  nspline1 = length(zero)-1,
  nspline2 = (length(zero)*2)-2)

#################
#Write the model#
#################

setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')
sink("combinedModel_grouped_spline.txt")
cat("
    model {
  
    #Model for the spline
    for (i in 1:n.Transect) {
      #spline effect
      abund.effect[i] <- eta[i]
      #random effect
      #abund.effect[i] <- randomS[i] + eta[i]
    }

    #random site
    for(i in 1:n.Transect){
      randomS[i] ~ dnorm(0,grid.tau)
    }
    grid.tau <- pow(grid.sd,-2)
    grid.sd ~ dunif(0,100)

    #the linear predictor
    eta <- X %*% b ## linear predictor for the spline
    
    ## Parametric effect priors
    b[1] ~ dnorm(1,0.001)
    
    ## prior for s(x,y)... 
    K1 <- S1[1:nspline1,1:nspline1] * lambda[1]  + S1[1:nspline1,nspline:nspline2] * lambda[2]
    b[2:nspline] ~ dmnorm(zero[2:nspline],K1) 
    
    ## smoothing parameter priors 
    for (i in 1:2) {
      rho[i] ~ dunif(-12,12)
      lambda[i] <- exp(rho[i])
    }
    
    #Line transect data:
    
    # MODEL GROUP SIZE
    intercept.groupsize ~ dnorm(0,0.001)

    #fixed year effect
    for(i in 1:2){
      yearEffect.gs[i] ~ dnorm(0,0.001)
    }
    yearEffect.gs[3] <- 0

    #random grid effect
    for(i in 1:n.Transect){
      gridEffect[i] ~ dnorm(0,grid.gs.tau)
    }
    grid.gs.tau <- pow(grid.gs.sd,-2)
    grid.gs.sd ~ dunif(0,100)
    
    for(i in 1:n.Detections){
      d.Groupsize[i] ~ dpois(exp.GroupSize[i])
      log(exp.GroupSize[i]) <- intercept.groupsize + 
                                  gridEffect[detectionSites[i]] +
                                  yearEffect.gs[d.Year[i]]
    }

    #get average per transect and year
    for(k in 1:n.TransectYrs){
      for(i in 1:n.Detections){
        grp.size[i,k] <-  exp.GroupSize[i] * TransYrIdx[i,k]
        }
      grp.jt[ty.combos[k,1], ty.combos[k,2]] <- sum(grp.size[,k])/max(1,sum(TransYrIdx[,k]))  
    }
    
    for(j in 1:n.Transect){
      for(t in 1:n.Yrs){
        grp[j,t] <- ifelse(equals(grp.jt[j,t],0),exp(intercept.groupsize+
                                                      yearEffect.gs[t]+
                                                      gridEffect[j]),grp.jt[j,t])
      }
    }

    # MODEL NUMBER OF GROUPS DETECTED
    
    # PRIORS

    #fit model
    for (j in 1:n.Transect){
      for (t in 1:n.Yrs){
        n[j,t] ~ dpois(nHat[j,t]) 
    
        #work out area of each cell that is surveyed
        surveyArea[j,t]<-((transectLengths[j,t]/1000)*(ESW[j,t]/1000)*2)/9
    
        #relate fraction seen to fraction available in whole cell
        nHat[j,t] <- expN[j,t]*surveyArea[j,t]
        log(expN[j,t]) <-  abund.effect[j] +yearEffect[t]
      }
    }
    
    #fixed year effect
    for(i in 1:2){
      yearEffect[i] ~ dnorm(0,0.01)
    }
    yearEffect[3] <- 0
    
    #predict density (number of indivdiuals per grid cell) for all sites 
    #using the abundance effect and average group size
    #and the intercept of the line transect model for number of groups
    for(j in 1:n.Transect){
        for(t in 1:n.Yrs){
          Density.jt[j,t] <- expN[j,t] * grp[j,t]
        }
      Density[j] <- mean(Density.jt[j,])
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
    b.d.0 ~ dunif(4,5)

    #random transect effect
    for(i in 1:n.d.transectId){
      random.transect[i] ~ dnorm(0,random.transect.tau)
    }
    random.transect.tau <- pow(random.transect.sd,-2)
    random.transect.sd ~ dunif(0,100)
    
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
    for(j in 1:n.Transect){
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

source('/Users/dianabowler/Documents/NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("random.transect","b.d.0","mean.esw",
            "random.transect.sd","intercept.year","yearEffect",
            "intercept.lt","yearEffect.gs",
            "avDensity","totDensity","Density")

#inits <- function(){list(sigma=runif(1,80,150))}

ni<-200000
nb<-50000
out1 <- jags(bugs.data, inits=NULL, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni,parallel=TRUE)

print(out1,2)
traceplot(out1)
save(out1,file="out1_final_ltonly_spline.RData")

########################################################################################

#Plotting predictions
#load("out1_final_ltonly.RData")
plot(bugs.data$n,out1$mean$nHat)
setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')

myGridDF3km$fits <- NA
myGridDF3km$fits[sites.lt]<-out1$mean$Density

#groupSizes[is.na(groupSizes)]<-0
#groupsizeMeans<-apply(groupSizes,1,mean)
#myGridDF3km$fits[sites.lt]<-groupsizeMeans

#totalsMeans<-apply(totalsInfo,1,mean)
#myGridDF3km$fits[sites.lt]<-totalsMeans

myGridDF3km$x <- myGridDF3km$x + medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y + medGridDF3km.y

predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
#png(file="combinedModel_grouped_spline_ltonly.png")
plot(predRaster)
plot(sws,add=T)
#dev.off()



#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline.RData")

#save similations for Bayesian p-value
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#########################################################################################

