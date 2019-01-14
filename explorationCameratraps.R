##################################################
#Plotting the camera trap data and basic analysis#
#################################################

source('C:/Users/diana.bowler/OneDrive - NINA/EldsDeer Population Assessment/GitFiles/formattingCameratraps.R')

#merge files of all deer observations and camera trap info data
deer<-merge(deer,locations,by=c("ID","Year"))
nrow(deer)#282

#get sum of deer at each lat and lon
deer<-ddply(deer,.(LATITUDE..N.,LONGITUDE..E.),summarise,numberDeer=sum(numberDeer))
coordinates(deer)<-deer[,c("LONGITUDE..E.","LATITUDE..N.")]
proj4string(deer)<-CRS("+proj=longlat +datum=WGS84")

#set locations at coordinates
coordinates(locations)<-locations[,c("LONGITUDE..E.","LATITUDE..N.")]
proj4string(locations)<-CRS("+proj=longlat +datum=WGS84")

#plotting
plot(sws)
plot(transects,add=T)
plot(locations,add=T,col="gray")
plot(deer,add=T,col="red")

##############################################################################################

######################################################################
#using basic GLMS#####################################################
#study effects of habitat on number of observations at each site/year#
######################################################################

deer<-merge(deer,locations,by=c("ID","Year"),all.y=T)

#get number of observations of deer at each camera
nuObs<-ddply(deer,.(ID),summarise,numberDeer=sum(numberDeer,na.rm=T))
nuObs<-merge(nuObs,locations,by="ID")

#get length of time each camera trap was active
nuObs$Days<-as.numeric(nuObs$TAKE.DOWN.DATE-nuObs$SETUP.DATE)

#some are missing because the take down date is missing
tapply(nuObs$Days[!is.na(nuObs$Days)],nuObs$Year[!is.na(nuObs$Days)],median)
2014 2016 
10   38

#assume they were active for 38 days
nuObs$TAKE.DOWN.DATE[is.na(nuObs$TAKE.DOWN.DATE)&nuObs$Year==2016]<-nuObs$SETUP.DATE[is.na(nuObs$TAKE.DOWN.DATE)&nuObs$Year==2016]+38
nuObs$TAKE.DOWN.DATE[is.na(nuObs$TAKE.DOWN.DATE)&nuObs$Year==2014]<-nuObs$SETUP.DATE[is.na(nuObs$TAKE.DOWN.DATE)&nuObs$Year==2014]+10

#calculate the day difference again
nuObs$Days<-as.numeric(nuObs$TAKE.DOWN.DATE-nuObs$SETUP.DATE)

#add month
nuObs$Month<-month(nuObs$SETUP.DATE)
nuObs$Month[which(nuObs$Month==12)]<-0

#get rid of missing data
nuObs<-subset(nuObs,WATER.POINT!="")
nuObs<-subset(nuObs,SALT.LICK!="")

#remove those whete the camera was not operational on take down
nuObs<-subset(nuObs,!CAMERA.OPERATIONAL.TAKE.DOWN%in%c("FALSE","NO","No(Lost)"))

#########################################
#test the effects of habitat variables##
########################################

summary(nuObs)

#are any of the variables correlated?
qplot(sqrt(ForestCover),sqrt(Topography),data=nuObs)
hist(sqrt(nuObs$ForestCover))
hist(sqrt(nuObs$Topography))

#test effects on number of deer seen
glm1<-glm(numberDeer~sqrt(ForestCover)+WATER.POINT+SALT.LICK+Month,offset=log(Days),family=quasipoisson,data=nuObs)
summary(glm1)
#water point seems to be most useful

#test effects on whether or not deer were present
nuObs$presenceDeer<-ifelse(nuObs$numberDeer>0,1,0)
glm1<-glm(presenceDeer~sqrt(ForestCover)+WATER.POINT,offset=logit(Days),family=quasibinomial,data=nuObs)
summary(glm1)
#forest cover most useful
glm1<-glm(presenceDeer~sqrt(ForestCover)+WATER.POINT,offset=Days,family=quasibinomial(link="cloglog"),data=nuObs)
summary(glm1)

#if present, what predicts numbers seen
nuObsPresence<-subset(nuObs,presenceDeer==1)
glm1<-glm(numberDeer-1~sqrt(ForestCover)+WATER.POINT+SALT.LICK+Month,offset=log(Days),family=quasipoisson,data=nuObsPresence)
summary(glm1)
#water point

#plotting the data
tapply(nuObsPresence$numberDeer,nuObsPresence$WATER.POINT,mean)
g1<-qplot(WATER.POINT,log(numberDeer),data=nuObs,geom="boxplot")
g2<-qplot(ForestCover,log(numberDeer),data=nuObs)#maybe should treat
g3<-qplot(SALT.LICK,log(numberDeer),data=nuObs,geom="boxplot")
g4<-qplot(factor(Month),log(numberDeer),data=nuObs,geom="boxplot")

library(gridExtra)
grid.arrange(g1,g2,g3,g4,ncol=2)

######################################################################################################