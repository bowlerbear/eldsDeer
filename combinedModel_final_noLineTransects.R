#################################################################################
source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/eldsDeer/formattingcombinedModel.R')

#Use the JAGAM function to get the BUGS code for the spline
#http://www.petrkeil.com/?p=2385
setwd('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment')
library(mgcv)
jags.ready <- jagam(Deer.ct~s(x, y,k=10), 
                    data=subset(myGridDF3km,!is.na(Deer.ct)), 
                    family="binomial", 
                    sp.prior="log.uniform",
                    file="jagam.ct.txt")

model1 <- gam(Deer.ct ~ s(x, y,k=2), 
                    data=subset(myGridDF3km,!is.na(Deer.ct)), 
                    family="binomial")

myGridDF3km$fits<-NA
myGridDF3km$fits[sites.ct]<-model1$fitted.values
qplot(x,y,data=myGridDF3km,colour=fits,size=5)
#smooth terms don't explain anything.... neither does forest cover for that matter...

#get the data bits we need from jags data
X = jags.ready$jags.data$X
S1 = jags.ready$jags.data$S1
zero = jags.ready$jags.data$zero

######################################
#Compile the remaining data for model#
######################################

bugs.data<-list(#camera trap data
                site.ct = sites.ct,
                y.ct = y,
                n.CameraTrapSites = dim(y)[1],
                n.reps = dim(y)[3],
                n.1kmGrids = dim(y)[2],
                #covariates
                Forest=as.numeric(scale(log(myGridDF3km$ForestCover+1)))[sites.ct],
                Fields=as.numeric(scale(log(myGridDF3km$Fields+1)))[sites.ct],
                Military=myGridDF3km$MilitaryF[sites.ct],
                Villages=as.numeric(scale(log(myGridDF3km$Villages+1)))[sites.ct],
                Water=waterMatrix,
                Month=monthMatrix,
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                Forest_All=as.numeric(scale(log(myGridDF3km$ForestCover+1))),
                Military_All=myGridDF3km$MilitaryF,
                Villages_All=as.numeric(scale(log(myGridDF3km$Villages+1))),
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
    
    #spline model
    eta <- X %*% b ## linear predictor
    
    ## Parametric effect priors CHECK tau=1/5.9^2 is appropriate!
    for (i in 1:1) { b[i] ~ dnorm(0,0.01) }
    
    ## prior for s(x,y)... 
    K1 <- S1[1:9,1:9] * lambda[1]  + S1[1:9,10:18] * lambda[2]
    b[2:10] ~ dmnorm(zero[2:10],K1) 
    
    ## smoothing parameter priors CHECK...
    for (i in 1:2) {
    rho[i] ~ dunif(-3,3)
    lambda[i] <- exp(rho[i])
    }
    
    #(1) Camera trap data
    
    #Priors
    theta.water ~ dnorm(0,0.001)
    alpha.month ~ dnorm(0,0.001)

    #intercepts for each model
    intercept.ct ~ dunif(0,1)
    mean.p.ct ~ dunif(0,1)
    int.theta ~ dunif(0,1) 

    #Model (multi-scale occupancy model using cloglog so we are esssentially modelling abundance!)
    for (i in 1:n.CameraTrapSites) { 
    z.ct[i] ~ dbern(psi[i])
    
    #(1)
    #cloglog(psi[i]) <- cloglog(intercept.ct) + beta.forest * Forest[i] + beta.military * Military[i] +
    #                   beta.villages * Villages[i] + beta.fields * Fields[i]

    #(2)
    #cloglog(psi[i]) <- cloglog(intercept.ct) + beta.forest * Forest[i] + beta.military * Military[i] +
    #                   beta.villages * Villages[i] 


    #(3)
    cloglog(psi[i]) <- eta[i]
    
    #(4)
    #cloglog(psi[i]) <- cloglog(intercept.ct) + beta.forest * Forest[i] + beta.military * Military[i]
    
    #(5)
    #cloglog(psi[i]) <- cloglog(intercept.ct) + beta.forest * Forest[i]


    #availbility for detection model
    for (j in 1:n.1kmGrids){
    
      a[i,j] ~ dbern(mu.a[i,j])
      mu.a[i,j] <- z.ct[i] * theta[i,j]
      cloglog(theta[i,j]) <- cloglog(int.theta) + theta.water * Water[i,j]
    
    #detection submodel 
    for (k in 1:n.reps) { # Loop over replicate surveys
      y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
      mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
      cloglog(p.ct[i,j,k]) <- cloglog(mean.p.ct) + alpha.month * Month[i,j,k]
    }
    }
    }

    #calculate the Bayesian p-value
    e <- 0.0001
    for (i in 1:n.CameraTrapSites) { 
      for (j in 1:n.1kmGrids){
        for (k in 1:n.reps) {
    # Fit assessments: Chi-squared test statistic and posterior predictive check
    chi2[i,j,k] <- pow((y.ct[i,j,k] - mu.ct[i,j,k]),2) / (sqrt(mu.ct[i,j,k])+e)         # obs.
    mu.ct.new[i,j,k] ~ dbern(mu.ct[i,j,k])      # Replicate (new) data set
    chi2.new[i,j,k] <- pow((mu.ct.new[i,j,k] - mu.ct[i,j,k]),2) / (sqrt(mu.ct[i,j,k])+e) # exp
        }
      }
    }

    # Add up discrepancy measures for entire data set
    for (i in 1:n.CameraTrapSites) { 
      for (j in 1:n.1kmGrids){
          fit1[i,j] <- sum(chi2[i,j,])                     
          fit.new1[i,j] <- sum(chi2.new[i,j,])            
      }
        fit2[i] <- sum(fit1[i,])                     
        fit.new2[i] <- sum(fit.new1[i,]) 

    }
    fit <- sum(fit2[])                     # Omnibus test statistic actual data
    fit.new <- sum(fit.new2[])             # Omnibus test statistic replicate data

    #predict occupancy across all sites
    #for(i in 1:n.sites){
    #  cloglog(allPsi[i]) <- cloglog(intercept.ct) + beta.forest * Forest_All[i] + 
    #                    beta.military * Military_All[i]
    #}

    }
    ",fill = TRUE)
sink()

###############
#Run the model#
###############

source('C:/Users/diana.bowler/OneDrive - NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("rho","intercept.ct","mean.p.ct","beta.forest","beta.military","beta.fields","beta.villages",
            "psi","theta.water","alpha.month","fit","fit.new","allPsi")

zst.ct <- apply(y, 1, max, na.rm=T)
ast <-apply(y,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

inits <- function(){list(z.ct = zst.ct,
                         a = ast)}

ni<-50000
nb<-10000
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
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,1))
plot(sws,add=T)
dev.off()

#pulling out all coeffients of interest
coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline_noLineTransect(3).RData")

#save similations
simslistDF<-list(fit=out1$sims.list$fit,fit.new=out1$sims.list$fit.new)
mean(simslistDF$fit.new>simslistDF$fit)
#1, 0.4838833
#2, 0.4879333
#3, 0.48675
#4, 0.4910333

########################################################################################

#for predictions across the whole area
myGridDF3km$fits<-out1$mean$allPsi
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
predRaster<-mask(predRaster,sws)
png(file="combinedModel_grouped_spline_noLineTransect(4).png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,1))
plot(sws,add=T)
dev.off()


coefTable<-data.frame(out1$summary)
coefTable$Parameter<-row.names(out1$summary)
save(coefTable,file="coefTable_combinedModel_grouped_spline_noLineTransect(4).RData")
#########################################################################################