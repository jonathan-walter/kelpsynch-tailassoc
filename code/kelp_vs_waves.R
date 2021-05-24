## Study the relationship between kelp loss and wave height

rm(list=ls())


library(wsyn)
library(ncf)
library(ecodist)
library(RColorBrewer)
library(nlme)
library(fields)
library(copula)
library(mvtnorm)
library(matrixcalc)
library(Matrix)
library(grDevices)
library(mgcv)

setwd("/Users/jonathanwalter/GitHub/kelp-synchrony/TailAssociation")

wavesq<-read.csv("wave_Quarterly_CleanedBasic.csv", header=F)
kelpq<-read.csv("kelp_Quarterly_CleanedBasic.csv", header=F)
no3<-as.matrix(read.csv("NO3_Annual_CleanedBasic.csv", header=F))
quarters<-read.csv("Quarters_CleanedBasic.csv")
locs<-read.csv("Locs_CleanedBasic.csv")

years<-unique(quarters$Year)

locs$Region<-rep(NA,nrow(locs))
locs$Region[1:242]<-1
locs$Region[243:nrow(locs)]<-2

quarters$waveyear<-c(rep(1986,2),rep(1987:2018,each=4),rep(2019,2))

waves.max<-matrix(NA, nrow(wavesq), length(years))
waves.mean<-matrix(NA, nrow(wavesq), length(years))
waves.Q4<-matrix(NA, nrow(wavesq), length(years))
kelp.loss<-matrix(NA, nrow(kelpq), length(years))

for(ii in 1:nrow(locs)){
  
  for(yy in 1:(length(years)-1)){
    
    if(any(is.na(kelpq[ii, quarters$Year==years[yy]]))){next}
    if(any(is.na(wavesq[ii, quarters$Year==years[yy]]))){next}
    
    waves.max[ii,yy]<-max(wavesq[ii, quarters$waveyear==years[yy]])
    waves.mean[ii,yy]<-mean(as.numeric(wavesq[ii, quarters$waveyear==years[yy]]))
    waves.Q4[ii,yy]<-wavesq[ii, quarters$waveyear==years[yy] & quarters$Quarter==4]
    kelp.loss[ii,yy]<-((kelpq[ii, quarters$waveyear==years[yy] & quarters$Quarter==3] 
                       - kelpq[ii, quarters$waveyear==years[yy] & quarters$Quarter==1])
                       / kelpq[ii, quarters$waveyear==years[yy] & quarters$Quarter==3])
    
  }
}

pdf("kelploss_vs_calmness.pdf", height=9, width=6.5)

par(mfrow=c(2,1))

plot(c(waves.max)*-1,c(kelp.loss), ylab="Proportional kelp loss", xlab="Minimum wave calmness"
     ,ylim=c(0,1), main="All sites")
plot(c(waves.mean)*-1, c(kelp.loss), ylab="Proportional kelp loss", xlab="Mean wave calmness"
     ,ylim=c(0,1), main="All sites")

plot(c(waves.max[locs$Region==1,])*-1,c(kelp.loss[locs$Region==1,]), ylab="Proportional kelp loss", xlab="Minimum wave calmness"
     ,ylim=c(0,1), main="Central coast")
plot(c(waves.mean[locs$Region==1,])*-1, c(kelp.loss[locs$Region==1,]), ylab="Proportional kelp loss", xlab="Mean wave calmness"
     ,ylim=c(0,1), main="Central coast")

plot(c(waves.max[locs$Region==2,])*-1,c(kelp.loss[locs$Region==2,]), ylab="Proportional kelp loss", xlab="Minimum wave calmness"
     ,ylim=c(0,1), main="Southern coast")
plot(c(waves.mean[locs$Region==2,])*-1, c(kelp.loss[locs$Region==2,]), ylab="Proportional kelp loss", xlab="Mean wave calmness"
     ,ylim=c(0,1), main="Southern coast")

dev.off()


## statistical modelling

kelp.Q2<-as.matrix(kelpq[, quarters$Quarter==2])
kelp.Q3<-as.matrix(kelpq[, quarters$Quarter==3])

moddat<-data.frame(kelp.loss=c(kelp.loss)
                   , min.calm=-1*c(waves.max)
                   , mean.calm=-1*c(waves.mean)
                   , kelp.Q2=c(kelp.Q2)
                   , kelp.Q3=c(kelp.Q3)
                   , mean.no3=c(no3)
                   , lat=locs$Lat
                   , lon=locs$Lon
                   )

moddat<-moddat[complete.cases(moddat),]
moddat<-moddat[!is.infinite(kelp.loss),]
moddat<-moddat[moddat$kelp.loss >= 0,]
#moddat$kelp.loss[moddat$kelp.loss < 0] <-0

mod<-gam(kelp.loss ~ s(mean.no3) + s(kelp.Q3) + s(mean.calm), data=moddat, gamma=1.4, family=gaussian(link=identity))
gam.check(mod)
summary(mod)

plot(mod, select=3)

pdf("kelploss_gam.pdf", width=3.25, height=8)
par(mfrow=c(3,1), mar=c(4.1,4.1,1.1,1.1))
plot(mod)
dev.off()