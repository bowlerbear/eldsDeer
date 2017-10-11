##############################################################
#This file formats the line transect data
#and then fits the model using the continuous distance data
#plus allowing dynanics over time
#this model is written in BUGS within the R script
#############################################################

#get the formatted line transect data file
source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formatLineTransect_forBUGS.R')
#this leads to the creation of "datafile"

#Notes:
#a few data points are missing from the year 2000 in the processed data sheet
#and site 24 in 1999  
#it appears these data sheets were lost

#ignore transects 17 onwards? these transects have hardly any observations
datafile<-subset(datafile,Transect<17)

#truncate the distance data
datafile$GroupSize[datafile$Distance>250]<-0
datafile$Obs[datafile$Distance>250]<-0
datafile$Distance[datafile$Distance>250]<-NA

#re-center the variables
datafile$Transect<-datafile$Transect-min(datafile$Transect)+1
datafile$T<-datafile$Year-min(datafile$Year)+1

#get info on transect lengths for each year and transect
transectDistances<-acast(datafile,Transect~T,value.var="TransectDist",fun=max)#there is only one value, just take the maximum
transectDistances[is.infinite(transectDistances)]<-1 #add in nominal distance when there was none

#add in empty data rows for transects 13 to 16 in 2001 when we dont have any data
datafileNewRow<-subset(datafile,Transect%in%c(13:16)&Year==2001)
datafileNewRow$Year<-2000
datafileNewRow$T<-3
datafileNewRow$TransectDist<-1
datafile<-rbind(datafile,datafileNewRow)

#re-ordering
datafile<-datafile[order(datafile$Year),]
datafile<-datafile[order(datafile$Transect),]
datafile<-datafile[order(datafile$Obs,decreasing=T),]

#Get number of groups seen per transect and each year
groupInfo<-acast(datafile,Transect~T,value.var="Obs",fun=sum,na.rm=T)

#get average group size per transect
groupSizes<-acast(datafile,Transect~T,value.var="GroupSize",fun=mean)
groupSizes[is.na(groupSizes)]<-0

#detectionInfo summary data
detectionInfo<-datafile[,c("Obs","Distance","GroupSize")]
detectionInfo<-subset(detectionInfo,Obs==1)

#look at frequency histogram of distribution distances
library(lattice)
histogram(detectionInfo$Distance)

#potential covariates for the process(P) and the detection (D) model
datafile$Region<-ifelse(datafile$Region=="Blue",0,1)
datafile$Observer<-paste(datafile$Year,datafile$Transect)
covarsP<-as.matrix(datafile[,c("Obs","T","Region","Transect","GroupSize")])
covarsD<-covarsP[covarsP[,1]==1,]

#TransectYears used in the model
transectInfo<-unique(datafile[,c("Transect","T")])
transectYears<-unique(transectInfo[,"T"])
transectTransect<-unique(transectInfo[,"Transect"])
transectRegion<-ifelse(transectTransect<7,0,1)

#index all data points by transect and year for use in the model
datafileObs<-subset(datafile,Obs==1)
TransYrIdx<-matrix(nrow=nrow(datafileObs),ncol=nrow(transectInfo))
TransYrIdx[]<-0
for(i in 1:nrow(datafileObs)){
  TransYrIdx[i,which(transectInfo$T==datafileObs$T[i]&
    datafile$Transect[i]==transectInfo$Transect)]<-1
}

#compile data
bugs.data<-list(nTransect = length(unique(datafile$Transect)), 
            nYrs = length(unique(datafile$T)),
            W = 250,
            L = transectDistances,
            n = groupInfo,
            groupSizes = groupSizes,
            N = nrow(detectionInfo),
            y = detectionInfo$Distance,
            gs = detectionInfo$GroupSize,
            covarsP = covarsP,
            covarsD = covarsD,
            #nTransectYrs = nrow(transectInfo),
            transectYears = transectYears,
            transectRegion = transectRegion,
            transectTransect = transectTransect,
            #ty.combos = transectInfo,
            #TransYrIdx = TransYrIdx,
            zeros.dist = rep(0,nrow(detectionInfo)))

#run model
cat("

model{
  
  # DETECTABILITY COMPONENT#####
  ##### PRIORS
  
  pi <- 3.141593

  # priors for fixed effect parms for half-normal detection parm sigma
  # intercept
  b.df.0 ~ dunif(0,20)        
  
  # fixed effects for region and group size
  b.df.GroupSize ~ dnorm(0,0.0001) #effect of group size on detection

  ##### Begin model for *all detections*
  
  for( i in 1:N){
    
    ##########
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    mu.df[i] <- b.df.0 + b.df.GroupSize*covarsD[i,5]

    ##########
    # estimate of sd and var, given coeffs above
    sig.df[i] <- exp(mu.df[i])
    sig2.df[i] <- sig.df[i]*sig.df[i]
    
    # effective strip width
    esw[i] <- sqrt(pi * sig2.df[i] / 2)
    f0[i] <- 1/esw[i] #assuming all detected on the line
    
    # LIKELIHOOD
    # using zeros trick
    y[i] ~ dunif(0,W)
    L.f0[i] <- exp(-y[i]*y[i] / (2*sig2.df[i])) * 1/esw[i] #y are the distances
    nlogL.f0[i] <-  -log(L.f0[i])
    zeros.dist[i] ~ dpois(nlogL.f0[i])
  }
 
  #Predict effective strip width for all sites (even without detection data)
  for (j in 1:nTransect){
    for (t in 1:nYrs){
    pred.sig[j,t] <- exp(b.df.0 + b.df.GroupSize*pred.gs[j,t])
    pred.sig2[j,t] <- pow(pred.sig[j,t],2)
    ESW.JT[j,t] <- sqrt(pi * pred.sig2[j,t] / 2)
    }
  }

  ######
  # MODEL GROUP SIZE FOR EACH GROUP DETECTED
  #####
  ##### PRIORS
  
    #intercept
    B.gs.0 ~ dnorm(0, 0.00001)
  
    #fixed effects
    B.gs.Region ~ dnorm(0,0.00001)
    B.gs.T ~ dnorm(0,0.00001) # time trend coefficient
    B.gs.Region.T ~ dnorm(0,0.00001) #interaction between time and region coefficient

    #random transect effects
    sd.gs.trans ~ dunif(0,10)
    tau.gs.trans <- pow(sd.gs.trans,-2)
    for(i in 1:nTransect){
      random.gs.trans[i] ~ dnorm(0,tau.gs.trans)
    }

    #random transect time slopes
    sd.gs.transtime ~ dunif(0,19)
    tau.gs.transtime <- pow(sd.gs.transtime,-2)
    for(j in 1:nTransect){
      random.gs.transtime[j] ~ dnorm(0,tau.gs.transtime)
    }

    #random time effects
    #sd.gs.time ~ dunif(0,19)
    #tau.gs.time <- pow(sd.gs.time,-2)
    #for(t in 1:nYrs){
    #  random.gs.time[t] ~ dnorm(0,tau.gs.time)
    #}

  ##### Begin model for all detections

  for (i in 1:N){
    #fit model
    mu.gs[i] <- exp(B.gs.0 + B.gs.T*covarsP[i,2] + B.gs.Region*covarsP[i,3] +
    B.gs.Region.T*covarsP[i,2]*covarsP[i,3] +random.gs.trans[covarsP[i,4]]+
    covarsP[i,2]*random.gs.transtime[covarsP[i,4]]) 
    gs[i] ~ dpois(mu.gs[i])
  }

  #get predicted group size for all years and transects
  for (j in 1:nTransect){
    for (t in 1:nYrs){
      pred.gs[j,t]<-exp(B.gs.0 + B.gs.T*transectYears[t] + B.gs.Region*transectRegion[j] +
      B.gs.Region.T*transectRegion[j]*transectYears[t]+random.gs.trans[j]+
      random.gs.transtime[j]*transectYears[t]) 
    }
  }

  ######
  # MODEL NUMBER OF GROUPS DETECTED
  #####
  ##### PRIORS
            
  #intercept        
  B.n.0 ~ dnorm(0,0.00001)
            
  #fixed effects        
  B.n.T ~ dnorm(0,0.00001) # time trend coefficient
  B.n.Region ~ dnorm(0,0.00001) # region coefficient
  B.n.Region.T ~ dnorm(0,0.00001) #interaction coefficient between time and region
    
  #random transect effects
  sd.trans ~ dunif(0,10)
  tau.trans <- pow(sd.trans,-2)
  for(j in 1:nTransect){
    random.trans[j] ~ dnorm(0,tau.trans)
  }

  #random effects
  sd.time ~ dunif(0,10)
    tau.time <- pow(sd.time,-2)
    for(t in 1:nYrs){
    random.time[t] ~ dnorm(0,tau.time)
    }

  #random transect slopes
    sd.transtime ~ dunif(0,10)
    tau.transtime <- pow(sd.trans,-2)
    for(j in 1:nTransect){
    random.transtime[j] ~ dnorm(0,tau.transtime)
    }

  for (j in 1:nTransect){
   for (t in 1:nYrs){
             
    # Poisson model to observed data   
    n[j,t] ~ dpois(nHat[j,t])
    EffectiveArea [j,t] <- (L[j,t]/1000) * (ESW.JT[j,t]/1000) * 2
    nHat[j,t] <- Density[j,t] * EffectiveArea [j,t]

    #density of groups           
    log(Density[j,t]) <- B.n.0 + B.n.T*transectYears[t] + 
              B.n.Region*transectRegion[j]+
              B.n.Region.T*transectYears[t]*transectRegion[j]+
              #random.time[t]+
              random.trans[j]+
              random.transtime[j]*transectYears[t]
   }
  }

  #Model overall density  
  for (t in 1:nYrs){  
    for (j in 1:nTransect){
      D.ty[j,t] <- (Density[j,t]*pred.gs[j,t])
    }
    D.tot[t] <- mean(D.ty[,t])
  }
  
  #aggregate over regions
  for (t in 1:nYrs){
          D.region[1,t] <- mean(D.ty[1:6,t])
          D.region[2,t] <- mean(D.ty[7:16,t])
  }

  #predictions just based on data
  for (t in 1:nYrs){  
    for (j in 1:nTransect){
      Ddata.ty[j,t] <- max(0,n[j,t]*groupSizes[j,t]*EffectiveArea [j,t])
    }
    Ddata.tot[t] <- mean(D.ty[,t])
  }
  
  #aggregate over regions
  for (t in 1:nYrs){
          Ddata.region[1,t] <- mean(Ddata.ty[1:6,t])
          Ddata.region[2,t] <- mean(Ddata.ty[7:16,t])
  }

  }
  ",fill=TRUE,file="continuous_Dynamic.txt")

  source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')
  
  params <- c("b.df.GroupSize","b.df.0","B.gs.0","B.n.0","B.gs.T","B.gs.Region","B.gs.Region.T","B.n.T",
              "B.n.Region","B.n.Region.T","Density","D.ty","D.tot","D.region","pred.gs",
              "Ddata.ty","Ddata.tot","Ddata.region")
  
  inits <- function(){list(b.df.0 = runif(1,2,5), 
                           B.gs.0 = runif(1,0.2,3),
                           B.n.0 = runif(1,0.5,5))}
  
  n.iter<-10000
  out1 <- jags(bugs.data, inits=inits, params, "continuous_Dynamic.txt", n.thin=nt,
               n.chains=nc, n.burnin=nb,n.iter=ni)
  
  print(out1,2)
  
  #to look at a certain parameter
  out1$summary[grepl("pred.gs",row.names(out1$summary)),]
  
  #look at histograms all of parameters
  myparams<-row.names(out1$summary)
  myparams<-sapply(myparams,function(x)strsplit(x,"\\[")[[1]][1])
  myData<-data.frame(Param=as.character(myparams),Coef=do.call(c,out1$mean))
  
  ggplot(myData)+
    geom_histogram(aes(x=Coef))+
    facet_wrap(~Param,scales="free")
  
  #############################################################################################
  
  ##########################
  #Plotting the predictions#
  ##########################
  setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
  #(1) get the predictions of densities per transect and time
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("D.ty",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:16,20)
  expectedDensities$Year<-rep(1:20,each=16)
  
  #plotting
  png("transect_ts.png")
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect,scales="free")
  dev.off()
  
  #(2) get the predictions of densities per time
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("D.tot",row.names(expectedDensities)),])
  expectedDensities$Year<-1998:2017
  
  #plotting
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+theme_bw()
  
  #(3) get the predictions of densities per region and time
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("D.region",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:2)
  expectedDensities$Year<-rep(1:20,each=2)
  
  #plotting
  png("regionts.png")
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect)
  dev.off()
  
  #(4) get the predictions of densities (number of groups) per time and transect
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("Density",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:16,20)
  expectedDensities$Year<-rep(1:20,each=16)
  
  #plotting
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect,scales="free")
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect)
  
  #(5) get the predictions of average group size per time and transect
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("pred.gs",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:16,20)
  expectedDensities$Year<-rep(1:20,each=16)
  
  #plotting
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect,scales="free")
  
  ggplot(expectedDensities)+geom_point(aes(x=Year,y=mean))+facet_wrap(~Transect,scales="free")
  
  #(6) plotting the predictions just based off the data - by region
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("Ddata.region",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:2)
  expectedDensities$Year<-rep(1:20,each=2)
  
  png("region_data_ts.png")
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect)
  dev.off()
  
  #(7) plotting the predictions just based off the data - by transect
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("Ddata.ty",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:16,20)
  expectedDensities$Year<-rep(1:20,each=16)
  
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+
    facet_wrap(~Transect,scales="free")
  
  ##############################################################################################
  
  #######################
  #Plotting the data#####
  #######################
  
  #number of groups
  nuGroups<-ddply(datafile,.(Year,Transect),summarise,nuGroups=sum(Obs))
  qplot(Year,nuGroups,data=nuGroups)+facet_wrap(~Transect)
  
  #average group size
  avGroups<-ddply(datafile,.(Year,Transect),summarise,meanGroupSize=mean(GroupSize))
  qplot(Year,meanGroupSize,data=avGroups)+facet_wrap(~Transect)
  
  #total individuals
  totalIndiv<-ddply(datafile,.(Year,Transect),summarise,total=sum(GroupSize))
  qplot(Year,total,data=totalIndiv)+facet_wrap(~Transect)
  
  ################################################################################################
  