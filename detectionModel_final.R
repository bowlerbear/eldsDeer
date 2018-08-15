#################################################################################

###################
#Retrive the data##
###################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')


####################################################################################################################

#plot the data to check for no mistakes

hist(detectionInfo$Distance)

####################################################################################################################

#####################################
#Compile the remaining data for model#
#####################################

bugs.data<-list( #line transect data
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
                d.GroupsizeS = as.numeric(scale(log(detectionInfo$GroupSize))),
                #d.Site = detectionInfo$Site,
                n.TransectYrs = nrow(transectInfo),
                ty.combos = transectInfo,
                TransYrIdx = TransYrIdx,
                zeros.dist = rep(0,nrow(detectionInfo)),
                transectLengths = mytransectLengths,
                n.transectId = length(unique(myGridDF3km$transectId[!is.na(myGridDF3km$Deer.lt)])),
                transectId = myGridDF3km$transectId[!is.na(myGridDF3km$Deer.lt)],
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)))

#################
#Write the model#
#################

setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
sink("detectionModel.txt")
cat("
    model {
    
    # DETECTABILITY COMPONENT#####

    pi <- 3.141593
    
    #Priors
    sigma ~ dunif(10,200)
    b.d.0 ~ dunif(0,10)
    b.group.size ~ dunif(0,1)

    #random transect effect
    for(i in 1:n.d.transectId){
        random.transect[i] ~ dnorm(0,random.transect.tau)
    }
    random.transect.tau <- pow(random.transect.sd,-2)
    random.transect.sd ~ dunif(0,10)

    ##### Begin model for *all detections*
    
    for( i in 1:n.Detections){
    
    #MODELS FOR HALF-NORMAL DETECTION PARAMETER SIGMA
    #mu.df[i] <- b.d.0
    mu.df[i] <- b.d.0 + random.transect[d.transectId[i]] + b.group.size * d.GroupsizeS[i]

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
    min.esw <- sqrt(pi * pow(exp(b.d.0),2) / 2) 
    for(j in 1:n.Transect){
      for(t in 1:n.Yrs){
        ESW[j,t] <- ifelse(equals(ESW.jt[j,t],0),min.esw,ESW.jt[j,t])
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

params <- c("b.group.size","b.d.0","random.transect.sd","sigma","min.esw","ESW")

ni<-10000
nb<-3000
out1 <- jags(bugs.data, inits=NULL, params, "detectionModel.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni,parallel=TRUE)

print(out1,2)
traceplot(out1)
#min.esw            120.78

########################################################################################
