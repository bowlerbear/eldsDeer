##############################################################
#This file formats the line transect data
#and then fits the model using the continuous distance data
#plus allowing dynanics over time
#this model is written in BUGS within the R script
#############################################################

#get the formatted line transect data file
source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/GitFiles/formatLineTransect_forBUGS.R')
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

#detectionInfo summary data
detectionInfo<-datafile[,c("Obs","Distance","GroupSize")]
detectionInfo<-subset(detectionInfo,Obs==1)

#look at frequency histogram of distribution distances
library(lattice)
histogram(detectionInfo$Distance)

#potential covariates for the process(P) and the detection (D) model
datafile$Observer<-paste(datafile$Year,datafile$Transect)
covarsP<-as.matrix(datafile[,c("Obs","T","Yellow","Transect","GroupSize")])
covarsD<-covarsP[covarsP[,1]==1,]

#TransectYears used in the model
transectInfo<-unique(datafile[,c("Transect","T")])
transectYears<-unique(transectInfo[,"T"])
transectTransect<-unique(transectInfo[,"Transect"])

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
            nTransectYrs = nrow(transectInfo),
            transectYears = transectYears,
            transectTransect = transectTransect,
            ty.combos = transectInfo,
            TransYrIdx = TransYrIdx,
            zeros.dist = rep(0,nrow(detectionInfo)))

#run model
cat("

model{
  
  # DETECTABILITY COMPONENT#####
  ##### PRIORS
  
  pi <- 3.141593
  a <- 0.147
  
  # priors for fixed effect parms for half-normal detection parm sigma
  # intercept
  b.df.0 ~ dunif(0,20)        
  
  # fixed effects for region and group size
  b.df.GroupSize ~ dnorm(0,0.0001) #effect of group size on detection
  #transect efffect 
  #b.df.Transect ~ dnorm(0,0.001)

  ##### Begin model for *all detections*
  
  for( i in 1:N){
    
    ##########
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    #mu.df[i] <- b.df.0 + b.df.GroupSize*covarsD[i,5]+b.df.Transect*covarsD[i,4]#no evidence for a transect effect after accouting for group size
    mu.df[i] <- b.df.0 + b.df.GroupSize*covarsD[i,5]

    ##########
    # estimate of sd and var, given coeffs above
    sig.df[i] <- exp(mu.df[i])
    sig2.df[i] <- sig.df[i]*sig.df[i]
    
    # effective strip width
    x[i] <- W/sqrt(2*sig2.df[i])   
    erf[i] <- sqrt(1-exp(-x[i]*x[i]*(4/pi + a*x[i]*x[i])/(1 + a*x[i]*x[i])))
    esw[i] <- sqrt(pi * sig2.df[i] / 2)  * erf[i]
    f0[i] <- 1/esw[i] #assuming all detected on the line
    
    # LIKELIHOOD
    # using zeros trick
    y[i] ~ dunif(0,W) 
    L.f0[i] <- exp(-y[i]*y[i] / (2*sig2.df[i])) * 1/esw[i] #y are the distances
    nlogL.f0[i] <-  -log(L.f0[i])
    zeros.dist[i] ~ dpois(nlogL.f0[i])
  }
 
  #get average esw per transect and year
  for(k in 1:nTransectYrs){
    for(i in 1:N){
    grp.ESW[i,k] <- esw[i] * TransYrIdx[i,k]
    }
    ESW[k] <- sum(grp.ESW[,k])  
    }
    
    #convert into an [j,t] array
    for(k in 1:nTransectYrs){
    ESW.JT[ty.combos[k,1], ty.combos[k,2]] <- ESW[k]/max(1,sum(TransYrIdx[,k]))
    }
    

  ######
  # MODEL GROUP SIZE FOR EACH GROUP DETECTED
  #####
  ##### PRIORS
  
  #intercept
  B.gs.0 ~ dnorm(0, 0.00001)
  
  #fixed effects
  B.gs.Transect ~ dnorm(0,0.00001)
  B.gs.T ~ dnorm(0,0.00001) # time trend coefficient
  B.gs.Transect.T ~ dnorm(0,0.00001) #interaction between time and region coefficient


    #random effects
    sd.gs.time ~ dunif(0,3)
    tau.gs.time <- pow(sd.gs.time,-2)
    for(t in 1:nYrs){
    random.gs.time[t] ~ dnorm(0,tau.gs.time)
    }

  ##### Begin model for all detections

  for (i in 1:N){
    
    # Transect + T
    mu.gs[i] <- exp(B.gs.0 + B.gs.T*covarsP[i,2] + B.gs.Transect*covarsP[i,4] +
    B.gs.Transect.T*covarsP[i,2]*covarsP[i,4]+random.gs.time[covarsP[i,2]]) 
    #mu.gs[i] <- exp(B.gs.0+random.time[covarsP[i,2]]) 
    gs[i] ~ dpois(mu.gs[i])
                      
  }

  #get average group size per transect and year
  for(k in 1:nTransectYrs){
    for(i in 1:N){
    grpSize[i,k] <- mu.gs[i] * TransYrIdx[i,k]
    }
    Es[k] <- sum(grpSize[,k])  
  }

  #convert into an [j,t] array
  for(k in 1:nTransectYrs){
    Es.gs[ty.combos[k,1], ty.combos[k,2]] <- Es[k]/max(1,sum(TransYrIdx[,k]))
  }

  ######
  # MODEL NUMBER OF GROUPS DETECTED
  #####
  ##### PRIORS
            
  #intercept        
  B.n.0 ~ dnorm(0,0.00001)
            
  #fixed effects        
  B.n.T ~ dnorm(0,0.00001) # time trend coefficient
  B.n.Transect ~ dnorm(0,0.00001) # region coefficient
  B.n.Transect.T ~ dnorm(0,0.00001) #interaction coefficient between time and region
    

  #random effects
  sd.time ~ dunif(0,3)
  tau.time <- pow(sd.time,-2)
  for(t in 1:nYrs){
    random.time[t] ~ dnorm(0,tau.time)
    }

  for (j in 1:nTransect){
   for (t in 1:nYrs){
             
    # Poisson model to observed data   
    n[j,t] ~ dpois(nHat[j,t])


    #density of groups           
    nHat[j,t] <- exp(B.n.0 + B.n.T*transectYears[t] + 
              B.n.Transect*transectTransect[j]+
              B.n.Transect.T*transectYears[t]*transectTransect[j]+
              random.time[t])

   }
  }

  #Model overall density  
  for (t in 1:nYrs){  
    for (j in 1:nTransect){
      D.ty[j,t] <- (Es.gs[j,t]*nHat[j,t])/(max(0.00001,(2*L[j,t]*ESW.JT[j,t])/1000))
      #Ddata[j,t] <- (groupSizes[j,t]*n[j,t])/(max(0.00001,(2*L[j,t]*ESW.JT[j,t])/1000))
    }
    Dtot[t] <- mean(D.ty[,t])
  }
  
  }
  ",fill=TRUE,file="continuous_Dynamic.txt")

  source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')
  
  params <- c("b.df.GroupSize","sig.df","b.df.0","B.gs.0","B.n.0","B.gs.T","B.gs.Transect","B.gs.Transect.T",
              "B.n.T","B.n.Transect","B.n.Transect.T","D.ty","Dtot")
  
  inits <- function(){list(b.df.0 = runif(1,2,5), 
                           B.gs.0 = runif(1,0.2,3),
                           B.n.0 = runif(1,0.5,5))}
  
  out1 <- jags(bugs.data, inits=inits, params, "continuous_Dynamic.txt", n.thin=nt,
               n.chains=nc, n.burnin=nb,n.iter=ni)
  
  print(out1,2)
  
  #sigma is about 80
  #average densities are
  
  #############################################################################################
  
  ##########################
  #Plotting the predictions#
  ##########################
  
  #(1) get the predictions of densities per transect and time
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("D.ty",row.names(expectedDensities)),])
  expectedDensities$Transect<-rep(1:16,18)
  expectedDensities$Year<-rep(1:18,each=16)
  
  #plotting
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+facet_wrap(~Transect,scales="free")
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.,colour=Transect))
  
  #(2) get the predictions of densities pertime
  expectedDensities<-out1$summary
  expectedDensities<-data.frame(expectedDensities[grepl("Dtot",row.names(expectedDensities)),])
  expectedDensities$Year<-1998:2015
  
  #plotting
  ggplot(expectedDensities)+geom_pointrange(aes(x=Year,y=mean,ymin=X2.5.,ymax=X97.5.))+theme_bw()
  #between 0.003 and 0.005
  
  ##############################################################################################
  