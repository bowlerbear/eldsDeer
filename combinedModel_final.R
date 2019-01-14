#################################################################################

###################
#Retrive the data##
###################

source('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')
old <- myGridDF3km
detectionInfo$Grid <- myGridDF3km$Grid[match(detectionInfo$Grid3km,myGridDF3km$Grid3km)]

####################################################################################################################

#plot the data to check for no mistakes

#plot all presences
myGrid3km[]<-0
myGrid3km[myGridDF3km$Grid3km]<-myGridDF3km$Deer
plot(myGrid3km)
plot(sws,add=T)

#plot line transect presences
myGrid3km[]<-0
myGrid3km[myGridDF3km$Grid3km]<-myGridDF3km$Deer.lt
plot(myGrid3km)

#plot camera trap presences
myGrid3km[]<-0
myGrid3km[myGridDF3km$Grid3km]<-myGridDF3km$Deer.ct
plot(myGrid3km)

#plot at 1km scale
myGrid3km[]<-rnorm(length(getValues(myGrid3km)))
myGrid[]<-0
myGrid[myGridDF$Grid1km]<-myGridDF$Deer.ct
plot(myGrid3km)
plot(myGrid,add=T)

#plot sampling regions
myGrid3km[]<-0
myGrid3km[myGridDF3km$Grid3km[myGridDF3km$Grid%in%sites.lt]]<-1
plot(myGrid3km)
myGrid3km[]<-0
myGrid3km[myGridDF3km$Grid3km[myGridDF3km$Grid%in%sites.ct]]<-1
plot(myGrid3km)

#plot total number of deer
myGridDF3km$DeerNu[myGridDF3km$Grid%in%sites.lt]<-apply(groupInfo,1,sum)
myGrid3km[]<-0
myGrid3km[myGridDF3km$Grid3km]<-myGridDF3km$DeerNu
plot(myGrid3km)
plot(sws,add=T)

####################################################################################################################

##############################################
#Sort out the code for the spline model#
##############################################
myGridDF3km <- old

medGridDF3km.x <- median(myGridDF3km$x)
medGridDF3km.y <- median(myGridDF3km$y)
myGridDF3km$x <- myGridDF3km$x - medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y - medGridDF3km.y

#Use the JAGAM function to get the BUGS code for the spline
setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')
library(mgcv)
df<-myGridDF3km 
m.gam<-gam(Deer ~ s(x, y,k=4),data=df, family=binomial)
gam.check(m.gam)
summary(m.gam)
#spline to residuals
df$resids <- residuals(m.gam)
m.gam<-gam(resids ~ s(x, y),data=df)
summary(m.gam)
######################################################################################################################

#for BUGS
jags.ready <- jagam(Deer ~ 1 + s(x, y,k=4), 
                    data=myGridDF3km, 
                    family=binomial, 
                    sp.prior="log.uniform",
                    file="jagam.txt")

#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero

########################################################################
#add NAs to the groupInfo data

new.my.n <-matrix(data=NA,nrow=nrow(myGridDF3km),ncol=3)
new.my.n[sites.lt,]<-groupInfo
new.transectlengths <-matrix(data=NA,nrow=nrow(myGridDF3km),ncol=3)
new.transectlengths[sites.lt,]<-mytransectLengths
new.transectlengths[is.na(new.transectlengths)]<-0
#new.transectlengths[new.transectlengths==0] <- 0.00000001

#index all data points
transectInfo2<-data.frame(expand.grid(Transect=1:nrow(myGridDF3km),T=1:3))
TransYrIdx2<-matrix(nrow=nrow(datafileObs),ncol=nrow(transectInfo2))
TransYrIdx2[]<-0
for(i in 1:nrow(datafileObs)){
  TransYrIdx2[i,which(datafileObs$T[i]==transectInfo2$T&
                       datafileObs$Transect[i]==transectInfo2$Transect)]<-1
}

############################################################################
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
  n = new.my.n,###groupInfo
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
  n.TransectYrs2 = nrow(transectInfo2),
  ty.combos = transectInfo,
  ty.combos2 = transectInfo2,
  TransYrIdx = TransYrIdx,
  TransYrIdx2 = TransYrIdx2,
  zeros.dist = rep(0,nrow(detectionInfo)),
  transectLengths = new.transectlengths,##mytransectLengths
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
  Month=monthMatrix-2,
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

    #Common model
    for (i in 1:n.sites) { #across all sites   
        
      #model with spline alone
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
  
      #Different observation models depending on the data:
  
      #(1) Camera trap data
      
      #Priors
      theta.lure ~ dnorm(0,0.01)
      theta.water ~ dnorm(0,0.01)
      alpha.month ~ dnorm(0,0.01)
      alpha.year ~ dnorm(0,0.01)      
      abund.slope ~ dnorm(0,0.01)
      alpha.forest ~ dnorm(0,0.01)
      alpha.field ~ dnorm(0,0.01)

      #intercepts for each model
      mean.p.ct ~ dunif(0,1)
      int.p <- cloglog(mean.p.ct)
      intercept.ct ~ dunif(0,1)
      int.ct <- cloglog(intercept.ct)
      int.theta ~ dunif(0,1) 
      int.t <- cloglog(int.theta)

    #Model (multi-scale occupancy model using cloglog so we are esssentially modelling abundance!)
    for (i in 1:n.CameraTrapSites) { 
        z.ct[i] ~ dbern(psi[i])
        cloglog(psi[i]) <- int.ct + abund.effect[site.ct[i]]
    
    #availbility for detection model
    for (j in 1:n.1kmGrids){
        a[i,j] ~ dbern(mu.a[i,j])
        mu.a[i,j] <- z.ct[i] * theta[i,j]
        #cloglog(theta[i,j]) <- int.t + theta.water * Water[i,j] + theta.lure * Lure[i,j]
        #cloglog(theta[i,j]) <- int.t + theta.water * Water[i,j]
        cloglog(theta[i,j]) <- int.t

    #detection submodel 
    for (k in 1:n.reps) { # Loop over replicate surveys
        y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
        mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
        #cloglog(p.ct[i,j,k]) <- int.p
        #cloglog(p.ct[i,j,k]) <- int.p + abund.slope * log(Density[site.ct[i]]) 
        cloglog(p.ct[i,j,k]) <- int.p + 
                                  #abund.slope * log(Density[site.ct[i]]) + 
                                  alpha.year * Year.mat[i,j,k]+ alpha.month * Month[i,j,k] 
                                  #alpha.forest * Forest[site.ct[i]] + alpha.field * Fields[site.ct[i]]
        }
      }
    }
    
    #(2) Line transect data:

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
    for(k in 1:n.TransectYrs2){
      for(i in 1:n.Detections){
        grp.size[i,k] <-  exp.GroupSize[i] * TransYrIdx2[i,k]
        }
      grp.jt[ty.combos2[k,1], ty.combos2[k,2]] <- sum(grp.size[,k])/max(1,sum(TransYrIdx2[,k]))  
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
          #Density.jt[j,t] <- exp(abund.effect[j]+yearEffect[t]) * grp[j,t]
          Density.jt[j,t] <- expN[j,t] * grp[j,t]
    }
          Density[j] <- mean(Density.jt[j,])
    }
    
    #get predicted total across whole area
    avDensity <- mean(Density)
    totDensity <- sum(Density)

    #calculate the Bayesian p-value
    e <- 0.0001
    #for(j in 1:n.Transect){
    #  for(t in 1:n.Yrs){
    #    # Fit assessments: Chi-squared test statistic and posterior predictive check
    #    chi2[j,t] <- pow((n[j,t] - nHat[j,t]),2) / (sqrt(nHat[j,t])+e)         # obs.
    #    nHat.new[j,t] ~ dpois(nHat[j,t])      # Replicate (new) data set
    #    chi2.new[j,t] <- pow((nHat.new[j,t] - nHat[j,t]),2) / (sqrt(nHat[j,t])+e) # exp
    #  }
    #}

    # Add up discrepancy measures for entire data set
    #for(t in 1:n.Yrs){
    #  fit.t[t] <- sum(chi2[,t])                     
    #  fit.new.t[t] <- sum(chi2.new[,t])             
    #}

    #fit <- sum(fit.t[])                     # Omnibus test statistic actual data
    #fit.new <- sum(fit.new.t[])             # Omnibus test statistic replicate data

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
                #b.group.size * d.GroupsizeS[i]+
                #b.forest * d.Forest[i]+
                #b.field * d.Field[i]

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
    for(k in 1:n.TransectYrs2){
      for(i in 1:n.Detections){
        grp.ESW[i,k] <- esw[i] * TransYrIdx2[i,k]
      }
        ESW.jt[ty.combos2[k,1], ty.combos2[k,2]] <- sum(grp.ESW[,k])/max(1,sum(TransYrIdx2[,k]))  
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

params <- c("b.d.0",
            "b.group.size","mean.esw","b.forest","b.field",
            "yearEffect.gs","yearEffect",
            "intercept.groupsize",
            "beta.forest","beta.military",
            "beta.villages","beta.fields",
            "abund.slope","alpha.forest","alpha.field",
            "alpha.month","theta.water","theta.lure","alpha.year",
            "gs.forest","gs.military","gs.field",
            "totDensity","Density")

zst.ct <- apply(bugs.data$y.ct, 1, max, na.rm=T)
ast <-apply(bugs.data$y.ct,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

inits <- function(){list(b.d.0 = runif(1,4.3,4.6),
                        z.ct = zst.ct, 
                         a = ast)}

ni<-200000
nb<-50000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni,parallel=TRUE)

print(out1,2)
traceplot(out1)

save(out1,file="out1_final_all_v4.RData")
save(out1,file="out1_final_all_sigs_final_inclYear2.RData")#with nas, 2
save(out1,file="out1_final_all_onlysigs_final_inclYear2.RData")#with nas, 2
save(out1,file="out1_final_splineonly_inclYear2.RData")#with nas, 2

########################################################################################

#Plotting predictions
myGridDF3km$x <- myGridDF3km$x + medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y + medGridDF3km.y

setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment")
load("out1_final_all_sigs_final_inclYear2.RData")
myGridDF3km$fits<-log(out1$mean$Density+1)
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
png(file="combinedModel.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,5.9))
plot(sws,add=T)
dev.off()

#fit for the spline
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment")
load("out1_final_splineonly_inclYear2.RData")
myGridDF3km$fits<-log(out1$mean$Density+1)
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
png(file="combinedModel_splineonly.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,5.9))
plot(sws,add=T)
dev.off()

#save similations for Bayesian p-value
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#Forest, 0.4820952
#ForestF, 0.4891429
#log F, 0.4644762

#########################################################################################