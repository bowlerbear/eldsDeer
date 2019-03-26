#################################################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

#Use the JAGAM function to get the BUGS code for the spline
#http://www.petrkeil.com/?p=2385
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
library(mgcv)
jags.ready <- jagam(Deer.ct~1+s(x, y), 
                    data=subset(myGridDF3km,!is.na(Deer.ct)), 
                    family="binomial", 
                    sp.prior="log.uniform",
                    file="jagam.ct.txt")

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
                #covariates
                Forest=as.numeric(scale(log(forestcover+1)))[sites.ct],
                ForestB=forestcoverB,
                Fields=as.numeric(scale(sqrt(fields+1))),
                Military=ifelse(military>0,1,0)[sites.ct],
                Villages=as.numeric(scale(sqrt(villages+1))),
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
    
    #(1) Camera trap data
    
    #Priors
    mean.p.ct ~ dunif(0, 1)
    alpha.ct <- cloglog(mean.p.ct)
    theta.water ~ dnorm(0,0.001)
    alpha.month ~ dnorm(0,0.001)
    intercept.ct ~ dnorm(0,0.001)

    # Intercepts availability probability
    for(j in 1:n.1kmGrids){
      int.thetaT[j] ~ dunif(0,1)
      int.theta[j] <- cloglog(int.thetaT[j])
    }
  
    #spline model
    eta <- X %*% b ## linear predictor

    ## prior for s(x,y)... 
    K1 <- S1[1:29,1:29] * lambda[1]  + S1[1:29,30:58] * lambda[2]
    b[1:29] ~ dmnorm(zero[1:29],K1) 
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
    rho[i] ~ dunif(-3,3)
    lambda[i] <- exp(rho[i])
    }

    #Model (multi-scale occupancy model using cloglog so we are esssentially modelling abundance!)
    for (i in 1:n.CameraTrapSites) { 
    z.ct[i] ~ dbern(psi[i])
    #(1)
    #logit(psi[i]) <- intercept.ct + beta.forest * Forest[i] 
    #(3)
    cloglog(psi[i]) <- intercept.ct + eta[i]
    
    #availbility for detection model
    for (j in 1:n.1kmGrids){
    a[i,j] ~ dbern(mu.a[i,j])
    mu.a[i,j] <- z.ct[i] * theta[i,j]
    cloglog(theta[i,j]) <- int.theta[j] + theta.water * Water[i,j]
    
    #detection submodel 
    for (k in 1:n.reps) { # Loop over replicate surveys
    y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
    mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
    cloglog(p.ct[i,j,k]) <- alpha.ct + alpha.month * Month[i,j,k]
    }
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

params <- c("intercept.ct","mean.p.ct","beta.forest","beta.military","psi",
            "theta.water","alpha.month")

zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

inits <- function(){list(z.ct = zst.ct,
                         a = ast)}

ni<-20000
nb<-5000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni)

print(out1,2)
traceplot(out1)

library(ggmcmc)   
ggs_traceplot(ggs(as.mcmc(out1)))

########################################################################################

#Plotting predictions

setwd("C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment")
myGridDF3km$fits<-NA
myGridDF3km$fits[sites.ct]<-out1$mean$psi
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_noLineTransect(3).png")
plot(predRaster)
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline_noLineTransect(3).RData")

#########################################################################################