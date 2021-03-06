#################################################################################

###################
#Retrive the data##
###################

source('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')
old <- myGridDF3km
####################################################################################################################
detectionInfo$Grid <- myGridDF3km$Grid[match(detectionInfo$Grid3km,myGridDF3km$Grid3km)]

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
df <- myGridDF3km
df$Deer.lt[is.na(df$Deer.lt)] <- df$Deer.ct[is.na(df$Deer.lt)]
m.gam<-gam(Deer.lt ~ s(x, y,k=4),data=df, family=poisson)
gam.check(m.gam)
summary(m.gam)
#spline to residuals
df$resids <- residuals(m.gam)
m.gam<-gam(resids ~ s(x, y),data=df)
summary(m.gam)

jags.ready <- jagam(Deer.lt~ 1 + s(x, y,k=4), 
                    data=df, 
                    family=poisson, 
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
#new.transectlengths[new.transectlengths==0] <- 0.00000001

#index all data points
transectInfo<-data.frame(expand.grid(Transect=1:nrow(myGridDF3km),T=1:3))
TransYrIdx<-matrix(nrow=nrow(datafileObs),ncol=nrow(transectInfo))
TransYrIdx[]<-0
for(i in 1:nrow(datafileObs)){
  TransYrIdx[i,which(datafileObs$T[i]==transectInfo$T&
                       datafileObs$Transect[i]==transectInfo$Transect)]<-1
}

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
  n = new.my.n,
  n.Detections = nrow(detectionInfo),
  n.detectionSites = length(unique(detectionInfo$Site)),
  detectionSites = detectionInfo$Transect,
  detectionGrids = detectionInfo$Grid,
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
  transectLengths = new.transectlengths,
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
    #Model of factors affecting abundance
    beta.forest ~ dnorm(0,0.01)  
    beta.military ~ dnorm(0,0.01)
    beta.fields ~ dnorm(0,0.01)
    beta.villages ~ dnorm(0,0.01)
    
    #random site
    for(i in 1:n.sites){
      randomS[i] ~ dnorm(0,grid.tau)
    }
    grid.tau <- pow(grid.sd,-2)
    grid.sd ~ dunif(0,10)

    #Common model
    for (i in 1:n.sites) { #across all sites   
    
    #model with the spline
      #abund.effect[i] <- eta[i] 

      #model with covariates and the spline
      #beta.fields * Fields[i]
      abund.effect[i] <- beta.forest * Forest[i] + beta.military * MilitaryN[i] +
                        eta[i]    
    }
    
    #Model for the spline
    #the linear predictor
    eta <- X %*% b ## linear predictor for the spline
    
    ## Parametric effect priors
    b[1] ~ dnorm(1,0.01)
    
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
    intercept.groupsize ~ dnorm(0,0.01)
    
    #effect of covariates
    gs.forest ~ dnorm(0,0.01)#no effect
    gs.military ~ dnorm(0,0.01)#no effect
    gs.field ~ dnorm(0,0.01)#no effect

    #fixed year effect
    for(i in 1:2){
      yearEffect.gs[i] ~ dnorm(0,0.01)
    }
    yearEffect.gs[3] <- 0
    
    #random grid effect
    for(i in 1:n.sites){
      gridEffect[i] ~ dnorm(0,grid.gs.tau)
    }
    grid.gs.tau <- pow(grid.gs.sd,-2)
    grid.gs.sd ~ dunif(0,10)
    
    for(i in 1:n.Detections){
      d.Groupsize[i] ~ dpois(exp.GroupSize[i])
      log(exp.GroupSize[i]) <- intercept.groupsize + 
                                  #gs.forest * d.Forest[i] + 
                                  #gs.military * d.Military[i] + 
                                  #gs.field * d.Field[i] +
                                  gridEffect[detectionGrids[i]] +
                                  yearEffect.gs[d.Year[i]]
    }

    #get average per transect and year
    for(k in 1:n.TransectYrs){
      for(i in 1:n.Detections){
        grp.size[i,k] <-  exp.GroupSize[i] * TransYrIdx[i,k]
        }
      grp.jt[ty.combos[k,1], ty.combos[k,2]] <- sum(grp.size[,k])/max(1,sum(TransYrIdx[,k]))  
    }
    
    for(j in 1:n.sites){
      for(t in 1:n.Yrs){
        grp[j,t] <- ifelse(equals(grp.jt[j,t],0),exp(intercept.groupsize+
                                                      yearEffect.gs[t]+
                                                      gridEffect[j]),grp.jt[j,t])
      }
    }

    # MODEL NUMBER OF GROUPS DETECTED
    
    # PRIORS
    #fit model
    for (j in 1:n.sites){
      for (t in 1:n.Yrs){
        n[j,t] ~ dpois(nHat[j,t]) 
    
        #work out area of each cell that is surveyed
        surveyArea[j,t]<-((transectLengths[j,t]/1000)*(ESW[j,t]/1000)*2)/9
    
        #relate fraction seen to fraction available in whole cell
        nHat[j,t] <- expN[j,t]*surveyArea[j,t]
        log(expN[j,t]) <- abund.effect[j] + yearEffect[t]
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
    for(j in 1:n.sites){
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
    b.group.size ~ dnorm(0,0.01)
    b.field ~ dnorm(0,0.01)
    b.forest ~ dnorm(0,0.01)
    b.military ~ dnorm(0,0.01)

    #random transect effect
    for(i in 1:n.d.transectId){
    random.transect[i] ~ dnorm(0,random.transect.tau)
    }
    random.transect.tau <- pow(random.transect.sd,-2)
    random.transect.sd ~ dunif(0,10)
    
    ##### Begin model for *all detections*
    
    for(i in 1:n.Detections){
    
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    #mu.df[i] <- b.d.0
    mu.df[i] <- b.d.0 + random.transect[d.transectId[i]] 
    #mu.df[i] <- b.d.0 + random.transect[d.transectId[i]] +
              # b.group.size * d.GroupsizeS[i]+
              #  b.forest * d.Forest[i]+
              # b.field * d.Field[i]

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

source('/Users/dianabowler/Documents/NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("random.transect","b.d.0",
            "b.group.size","b.forest","b.field","b.military","mean.esw",
            "intercept.year","yearEffect",
            "intercept.lt","intercept.groupsize",
            "averagePa","beta.forest","beta.forest2","beta.military",
            "abund.slope","beta.year",
            "beta.villages","beta.fields",
            "alpha.month","theta.water","theta.lure",
            "gs.forest","gs.military","gs.field","yearEffect.gs",
            "totDensity","Density","grp","fit","fit.new")

#inits <- function(){list(sigma=runif(1,80,150))}

ni<-200000
nb<-50000
out1 <- jags(bugs.data, inits=NULL, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni,parallel=TRUE)

print(out1,2)
traceplot(out1)
save(out1,file="out1_final_ltonly_allVars.RData")
save(out1,file="out1_final_ltonly_sig.RData")
save(out1,file="out1_final_ltonly_onlysig.RData")
save(out1,file="out1_final_ltonly_spline_inclYear2.RData")

########################################################################################

myGridDF3km$x <- myGridDF3km$x + medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y + medGridDF3km.y

#Plotting predictions
setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')
load("out1_final_ltonly_onlysig.RData")
myGridDF3km$fits<-log(out1$mean$Density+1)
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
png(file="combinedModel_lt.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,5.9))
plot(sws,add=T)
dev.off()

#spline only
setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')
load("out1_final_ltonly_spline_inclYear2.RData")
myGridDF3km$fits<-NA
myGridDF3km$fits[sites.lt]<-log(out1$mean$Density[sites.lt]+1)
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
png(file="combinedModel_lt_spline.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,5.9))
plot(sws,add=T)
dev.off()

#save similations for Bayesian p-value
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#0.3349524, ForestF
#0.4086667, log Forest
#0.3629524, Forest
#########################################################################################


