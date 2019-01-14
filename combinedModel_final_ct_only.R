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

medGridDF3km.x <- median(myGridDF3km$x)
medGridDF3km.y <- median(myGridDF3km$y)
myGridDF3km$x <- myGridDF3km$x - medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y - medGridDF3km.y

#Use the JAGAM function to get the BUGS code for the spline
setwd('/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment')
library(mgcv)
df <- subset(myGridDF3km,!is.na(Deer.ct))
m.gam<-gam(Deer.ct ~ s(x, y,k=4),data=df, family=binomial)
gam.check(m.gam)
summary(m.gam)
#spline to residuals
df$resids <- residuals(m.gam)
m.gam<-gam(resids ~ s(x, y),data=df)
summary(m.gam)
######################################################################################################################

#for BUGS
jags.ready <- jagam(Deer.ct ~ 1 + s(x, y,k=4), 
                    data=subset(myGridDF3km,!is.na(Deer.ct)), 
                    family=binomial, 
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
                #total number of sites  
                n.sites = length(unique(myGridDF3km$Grid3km)),
                #covariates
                Forest=as.numeric(scale(sqrt(myGridDF3km$ForestCover+0.01))),
                Forest2=as.numeric(scale(myGridDF3km$ForestCover^2)),
                Fields=as.numeric(scale(sqrt(myGridDF3km$Fields+0.01))),
                MilitaryN=as.numeric(scale(sqrt(myGridDF3km$Military+0.01))),
                Villages=as.numeric(sqrt(myGridDF3km$Village)),
                Forest.ct=as.numeric(scale(sqrt(myGridDF3km$ForestCover+0.01)))[sites.ct],
                Fields.ct=as.numeric(scale(sqrt(myGridDF3km$Fields+0.01)))[sites.ct],
                MilitaryN.ct=as.numeric(scale(sqrt(myGridDF3km$Military+0.01)))[sites.ct],
                #covariates for trap locations
                Water=waterMatrix,
                Year.mat=yearMatrix,
                Lure=lureMatrix,
                Month=monthMatrix-2,
                Forest.mat=forestMatrix,
                Field.mat=fieldMatrix,
                #spatial covariates
                X = X,
                S1 = S1,
                zero = zero,
                nspline = length(zero),
                nspline1 = length(zero)+1,
                nspline2 = length(zero)*2)

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

    for (i in 1:n.CameraTrapSites) {   
      
      #model with spline alone
      #abund.effect[i] <- eta[i]

      #beta.fields * Fields.ct[i]
      #model with covariates and the spline
      abund.effect[i] <- beta.forest * Forest.ct[i] + beta.military * MilitaryN.ct[i] +
                        eta[i] 
       
    }
 
      #Model for the spline
      #the linear predictor
      eta <- X %*% b ## linear predictor for the spline
  
      ## prior for s(x,y)... 
      K1 <- S1[1:nspline,1:nspline] * lambda[1]  + S1[1:nspline,nspline1:nspline2] * lambda[2]
      b[1:nspline] ~ dmnorm(zero[1:nspline],K1)

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
      alpha.forest ~ dnorm(0,0.01)
      alpha.field ~ dnorm(0,0.01)
      yearEffect ~ dnorm(0,0.01)
  
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
        cloglog(psi[i]) <- int.ct + abund.effect[i] 
    
    #availbility for detection model
    for (j in 1:n.1kmGrids){
        a[i,j] ~ dbern(mu.a[i,j])
        mu.a[i,j] <- z.ct[i] * theta[i,j]
        cloglog(theta[i,j]) <- int.t
        #cloglog(theta[i,j]) <- int.t + theta.water * Water[i,j] + theta.lure * Lure[i,j]
        #cloglog(theta[i,j]) <- int.t + theta.water * Water[i,j]
        #cloglog(theta[i,j]) <- int.t + theta.water * Water[i,j] + theta.lure * Lure[i,j] +
        #                        alpha.forest * Forest.mat[i,j] + alpha.field * Field.mat[i,j]

    #detection submodel 
    for (k in 1:n.reps) { # Loop over replicate surveys
        y.ct[i,j,k] ~ dbern(mu.ct[i,j,k])    
        mu.ct[i,j,k] <- a[i,j] * p.ct[i,j,k]
        cloglog(p.ct[i,j,k]) <- int.p + alpha.month * Month[i,j,k] + alpha.year * Year.mat[i,j,k]
        #cloglog(p.ct[i,j,k]) <- int.p + alpha.month * Month[i,j,k] + alpha.year * Year.mat[i,j,k]+
                                #alpha.field * Fields.ct[i] + alpha.forest * Forest.ct[i] 
        }
      }
    }
    
    #predict across all sites
    for(i in 1:n.sites){
      cloglog(Density[i]) <-  int.ct + beta.forest * Forest[i] + 
                              beta.military * MilitaryN[i] 
    }


    }
    ",fill = TRUE)
sink()

###############
#Run the model#
###############

source('/Users/dianabowler/Documents/NINA/methods/models/bugsFunctions.R')

#https://www.r-bloggers.com/spatial-autocorrelation-of-errors-in-jags/

params <- c("beta.forest","beta.forest2","beta.military",
            "beta.villages","beta.fields","mean.p.ct",
            "alpha.month","theta.water","theta.lure","alpha.year","alpha.forest","alpha.field",
            "Density","psi","z.ct")

zst.ct <- apply(bugs.data$y.ct, 1, max, na.rm=T)
ast <-apply(bugs.data$y.ct,c(1,2),max,na.rm=T)
ast[is.infinite(ast)]<-0

inits <- function(){list(z.ct = zst.ct, 
                         a = ast,
                         mean.p.ct = 0.1)}

ni<-600000
nb<-200000
out1 <- jags(bugs.data, inits=inits, params, "combinedModel_grouped_spline.txt", n.thin=nt,
             n.chains=nc, n.burnin=nb,n.iter=ni,parallel=TRUE)

print(out1,2)
traceplot(out1)

save(out1,file="out1_final_all_v2.RData")#done
save(out1,file="out1_final_all_sigs_pre.RData")#done
save(out1,file="out1_final_all_sigs.RData")#done
save(out1,file="out1_final_splineonly.RData")#done

########################################################################################

setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment")
myGridDF3km$x <- myGridDF3km$x + medGridDF3km.x
myGridDF3km$y <- myGridDF3km$y + medGridDF3km.y

#Plotting predictions
load("out1_final_all_sigs.RData")
myGridDF3km$fits<-out1$mean$Density
myGridDF3km$fits[sites.ct]<-out1$mean$psi
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
png(file="combinedModel_ct.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,1))
plot(sws,add=T)
dev.off()

#fit for the spline
setwd("/Users/dianabowler/Documents/NINA/EldsDeer Population Assessment")
load("out1_final_splineonly.RData")
myGridDF3km$fits<-NA
myGridDF3km$fits[sites.ct]<-out1$mean$psi
predRaster<-rasterFromXYZ(myGridDF3km[,c("x","y","fits")])
png(file="combinedModel_ct_splineonly.png")
plot(predRaster,col=rev(topo.colors(50)),zlim=c(0,1))
plot(sws,add=T)
dev.off()

#########################################################################################