rm(list=ls())

library(wsyn)
library(ncf)
library(ecodist)
library(RColorBrewer)
library(nlme)
library(fields)
library(copula)
library(mvtnorm)
#library(matrixcalc)
library(Matrix)
library(grDevices)
library(rgdal)

setwd("~/GitHub/kelpsynch-tailassoc")

source("./code/partialSpearman.R")
source("./code/ncsurrog.R")

kelp.raw<-as.matrix(read.csv("./data/Kelp_Annual_CleanedBasic.csv", header=F))
no3.raw<-as.matrix(read.csv("./data/NO3_Annual_CleanedBasic.csv", header=F))
waves.raw<-as.matrix(read.csv("./data/wave_Annual_CleanedBasic.csv",header=F))
clims.raw<-read.csv("./data/Climind_Annual_CleanedBasic.csv")
quarters<-read.csv("./data/Quarters_CleanedBasic.csv")
locs<-read.csv("./data/Locs_CleanedBasic.csv")

#detrend the data
kelp<-kelp.raw
no3<-no3.raw
waves<-waves.raw
clims<-clims.raw

years<-unique(quarters$Year)

for(ii in 1:nrow(kelp)){
  
  kelp[ii,]<-residuals(lm(kelp.raw[ii,] ~ years))
  no3[ii,]<-residuals(lm(no3.raw[ii,] ~years))
  waves[ii,]<-residuals(lm(waves[ii,]~years))

}

clims$NPGO<-residuals(lm(clims$NPGO ~ years))

locs$Region<-rep(NA,nrow(locs))
locs$Region[1:242]<-1
locs$Region[243:nrow(locs)]<-2

# pdf("locs_by_cluster.pdf")
# plot(locs$Lon, locs$Lat, col=locs$Region, pch=19, cex=0.2)
# legend("topright",pch=19,col=1:2,legend=c("Central","Southern"))
# dev.off()


lb<-c(0,0.5)
ub<-c(0.5,1)

cormat.all<-cor(t(kelp), method="spearman")
filt.pos<-ifelse(cormat.all>=0,1,NA)

cormat.lb<-matrix(0,nrow(kelp),nrow(kelp))
cormat.ub<-matrix(0,nrow(kelp),nrow(kelp))
for(ii in 2:nrow(kelp)){
  for(jj in 1:(ii-1)){
    
    cormat.lb[ii,jj]<-partialSpearman(kelp[ii,],kelp[jj,],lb)
    cormat.ub[ii,jj]<-partialSpearman(kelp[ii,],kelp[jj,],ub)
    
  }
}

cormat.lb<-cormat.lb + t(cormat.lb)
cormat.ub<-cormat.ub + t(cormat.ub)

cormat.lb.filt<-cormat.lb*filt.pos
cormat.ub.filt<-cormat.ub*filt.pos

#distance decay --  does it depend on tail association?

splineFit<-function(distmat,zmat,nresamp=1000,quantiles=c(0,0.01,0.025,0.05,0.1,0.5,0.9,0.95,0.975,0.99,1)){
  triang<-lower.tri(distmat)
  distmat<-distmat
  xemp<-distmat[triang]
  yemp<-zmat[triang]
  drop.NaNs<-!is.na(yemp)
  dfs=sqrt(nrow(distmat))
  out<-list()
  
  emp.spline<-smooth.spline(xemp[drop.NaNs],yemp[drop.NaNs],df=dfs)
  out$emp.spline<-emp.spline
  
  resamp.splines<-matrix(NA, nrow=nresamp, ncol=length(emp.spline$y))
  for(ii in 1:nresamp){
    shuffle<-sample(1:nrow(distmat), size=nrow(distmat), replace=TRUE)
    xres<-distmat[shuffle,shuffle][triang]
    yres<-zmat[shuffle,shuffle][triang]
    drop.NaNs<-!is.na(yres)
    xres<-xres[drop.NaNs]
    yres<-yres[drop.NaNs]
    yres<-yres[!(xres==0)]
    xres<-xres[!(xres==0)]
    res.spline<-smooth.spline(xres,yres,df=dfs)
    resamp.splines[ii,]<-predict(res.spline, x=emp.spline$x)$y
  }
  out$resamp.splines<-resamp.splines
  out$spline.quantiles<-apply(resamp.splines,2,quantile,probs=quantiles)
  return(out)
}

geog.dist<-gcdist(locs$Lon,locs$Lat)

sncf.lb<-splineFit(geog.dist, cormat.lb.filt, quantiles=c(0.025,0.5,0.975))
sncf.ub<-splineFit(geog.dist, cormat.ub.filt, quantiles=c(0.025,0.5,0.975))


# plot(geog.dist[lower.tri(geog.dist)], cormat.lb[lower.tri(cormat.lb)], main="Lower tail",
#      xlab="Distance (Km)", ylab="Partial Spearman correlation", pch=20, cex=0.2, ylim=c(-0.5,0.5))
# lines(sncf.lb$emp.spline$x, sncf.lb$emp.spline$y, col="red")
# lines(sncf.lb$emp.spline$x,sncf.lb$spline.quantiles[1,], col="red", lty=2)
# #lines(sncf.lb$emp.spline$x,sncf.lb$spline.quantiles[2,], col="red", lty=2)
# lines(sncf.lb$emp.spline$x,sncf.lb$spline.quantiles[3,], col="red", lty=2)
# 
# plot(geog.dist[lower.tri(geog.dist)], cormat.lb[lower.tri(cormat.ub)], main="Upper tail", 
#      xlab="Distance (Km)", ylab="Partial Spearman correlation", pch=20, cex=0.2, ylim=c(-0.5,0.5))
# lines(sncf.ub$emp.spline$x, sncf.ub$emp.spline$y, col="blue")
# lines(sncf.ub$emp.spline$x,sncf.ub$spline.quantiles[1,], col="blue", lty=3)
# #lines(sncf.lb$emp.spline$x,sncf.ub$spline.quantiles[2,], col="blue", lty=2)
# lines(sncf.ub$emp.spline$x,sncf.ub$spline.quantiles[3,], col="blue", lty=3)

# pdf("dist_decay_allsites.pdf")
# plot(sncf.lb$emp.spline$x, sncf.lb$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
#      xlab="Distance (km)", ylab="Partial Spearman correlation", main="All coastline segments")
# lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[1,], col="red", lty=3)
# lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[3,], col="red", lty=3)
# lines(sncf.ub$emp.spline$x, sncf.ub$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[1,], col="blue", lty=3)
# lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[3,], col="blue", lty=3)
# legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
# dev.off()



#look at regions separately -------------------------------------------------------------------------------------------

geog.dist.c1<-geog.dist[locs$Region==1,locs$Region==1]
geog.dist.c2<-geog.dist[locs$Region==2,locs$Region==2]

cormat.lb.c1<-cormat.lb.filt[locs$Region==1,locs$Region==1]
cormat.ub.c1<-cormat.ub.filt[locs$Region==1,locs$Region==1]
cormat.lb.c2<-cormat.lb.filt[locs$Region==2,locs$Region==2]
cormat.ub.c2<-cormat.ub.filt[locs$Region==2,locs$Region==2]

sncf.lb.c1<-splineFit(geog.dist[locs$Region==1,locs$Region==1], cormat.lb.c1, quantiles=c(0.025,0.5,0.975))
sncf.ub.c1<-splineFit(geog.dist[locs$Region==1,locs$Region==1], cormat.ub.c1, quantiles=c(0.025,0.5,0.975))
sncf.lb.c2<-splineFit(geog.dist[locs$Region==2,locs$Region==2], cormat.lb.c2, quantiles=c(0.025,0.5,0.975))
sncf.ub.c2<-splineFit(geog.dist[locs$Region==2,locs$Region==2], cormat.ub.c2, quantiles=c(0.025,0.5,0.975))


# plot(geog.dist.c1[lower.tri(geog.dist)], cormat.lb.c1[lower.tri(cormat.lb)], main="C1 Lower tail",
#      xlab="Distance (Km)", ylab="Partial Spearman correlation", pch=20, cex=0.2, ylim=c(-0.5,0.5))
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$emp.spline$y, col="red")
# lines(sncf.lb.c1$emp.spline$x,sncf.lb.c1$spline.quantiles[1,], col="red", lty=2)
# #lines(sncf.lb$emp.spline$x,sncf.lb$spline.quantiles[2,], col="red", lty=2)
# lines(sncf.lb.c1$emp.spline$x,sncf.lb.c1$spline.quantiles[3,], col="red", lty=2)
# 
# plot(geog.dist.c1[lower.tri(geog.dist)], cormat.ub.c1[lower.tri(cormat.ub)], main="C1 Upper tail", 
#      xlab="Distance (Km)", ylab="Partial Spearman correlation", pch=20, cex=0.2, ylim=c(-0.5,0.5))
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$emp.spline$y, col="blue")
# lines(sncf.ub.c1$emp.spline$x,sncf.ub.c1$spline.quantiles[1,], col="blue", lty=3)
# #lines(sncf.lb$emp.spline$x,sncf.ub$spline.quantiles[2,], col="blue", lty=2)
# lines(sncf.ub.c1$emp.spline$x,sncf.ub.c1$spline.quantiles[3,], col="blue", lty=3)

# pdf("dist_decay_central.pdf")
# plot(sncf.lb.c1$emp.spline$x, sncf.lb.c1$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
#      xlab="Distance (km)", ylab="Partial Spearman correlation", main="Central California (c1)")
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[1,], col="red", lty=3)
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[3,], col="red", lty=3)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[1,], col="blue", lty=3)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[3,], col="blue", lty=3)
# legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
# dev.off()

# pdf("dist_decay_southern.pdf")
# plot(sncf.lb.c2$emp.spline$x, sncf.lb.c2$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
#      xlab="Distance (km)", ylab="Partial Spearman correlation", main="Southern California (c2)")
# lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[1,], col="red", lty=3)
# lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[3,], col="red", lty=3)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[1,], col="blue", lty=3)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[3,], col="blue", lty=3)
# legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
# dev.off()


## visualize differences between upper and lower tails---------------------------------------------

cormat.diff<-cormat.ub.filt-cormat.lb.filt

# pdf("hist_kelpsynch_taildiff.pdf")
# hist(cormat.diff[lower.tri(cormat.diff)])
# abline(v=median(cormat.diff[lower.tri(cormat.diff)], na.rm=T), lty=2, col="red")
# dev.off()

mean(cormat.diff[lower.tri(cormat.diff)], na.rm=T)

cormat.diff.nodiag<-cormat.diff
diag(cormat.diff.nodiag)<-NA

# pdf("synchmat_kelp_taildiff.pdf")
# image.plot(cormat.diff.nodiag)
# dev.off()



## test whether geographies of kelp synchrony differ between upper and lower tails ----------------
# mantel(as.dist(cormat.lb)~as.dist(cormat.ub))
# 
# set.seed(17)
# 
# ncsurr<-ncsurrog(m=t(kelp), numsurrog=1000, plotcheckloc=paste0(getwd(),"/plotcheck"))
# 
# empDiff<-sum(cormat.ub.filt-cormat.lb.filt, na.rm=T)
# empCorr<-cor(cormat.lb.filt[lower.tri(cormat.lb.filt)]
#              ,cormat.ub.filt[lower.tri(cormat.ub.filt)]
#              ,use="pairwise.complete.obs")
# 
# surrDiff<-rep(NA, dim(ncsurr)[3])
# surrCorr<-rep(NA, dim(ncsurr)[3])
# for(rep in 1:dim(ncsurr)[3]){
# 
#   N<-dim(ncsurr)[2]
# 
#   tmpdat<-t(ncsurr[,,rep])
#   tmpcormat<-cor(t(tmpdat), method="spearman")
#   tmpfilt<-ifelse(tmpcormat>=0,1,NA)
#   tmpcormat.lb<-matrix(0,N,N)
#   tmpcormat.ub<-matrix(0,N,N)
#   for(ii in 2:nrow(kelp)){
#     for(jj in 1:(ii-1)){
#       tmpcormat.lb[ii,jj]<-partialSpearman(tmpdat[ii,],tmpdat[jj,],lb)
#       tmpcormat.ub[ii,jj]<-partialSpearman(tmpdat[ii,],tmpdat[jj,],ub)
#     }
#   }
# 
#   tmpcormat.lb<-tmpcormat.lb + t(tmpcormat.lb) * tmpfilt
#   tmpcormat.ub<-tmpcormat.ub + t(tmpcormat.ub) * tmpfilt
# 
#   surrDiff[rep]<-sum(tmpcormat.ub-tmpcormat.lb, na.rm=T)
#   surrCorr[rep]<-1-cor(tmpcormat.lb[lower.tri(tmpcormat.lb)]
#                      ,tmpcormat.ub[lower.tri(tmpcormat.ub)]
#                      ,use="pairwise.complete.obs")
# }
# 
# rank(c(empDiff,surrDiff))[1]/1001 ## overall degree of tail dependence not significant
# rank(c(empCorr,surrCorr))[1]/1001 ## but geography significantly differs


## Measure average synchrony within a distance threshold

avg_by_dist<-function(cormat, dmat, thresh){
  
  diag(cormat)<-NA
  
  nbhd<-ifelse(dmat <= thresh, 1, 0)
  nbhd.cor<-cormat * nbhd
  return(apply(nbhd.cor, 1, mean, na.rm=T))
  
}

geog.dist<-ncf::gcdist(locs$Lon, locs$Lat)

coravg.lb.25km<-avg_by_dist(cormat.lb.filt, geog.dist, 25)
coravg.ub.25km<-avg_by_dist(cormat.ub.filt, geog.dist, 25)
coravg.diff.25km<-avg_by_dist(cormat.diff, geog.dist, 25)

#coravg.diff.25km<-abs(coravg.lb.25km)-abs(coravg.ub.25km)

nclasses=7
pal<-brewer.pal(nclasses,"RdYlBu")
pal2<-brewer.pal(nclasses, "Greens")



class_eq_int<-function(x, n){
  
  breaks<-seq(from=min(x, na.rm=T), to=max(x, na.rm=T)+1/length(x), length.out=n+1)
  class<-rep(NA, length(x))
  
  for(ii in 1:n){
    class[x>=breaks[ii] & x<breaks[ii+1]]<-ii
  }
  return(list(class=class, breaks=breaks))
}


class_eq_int_sym<-function(x, n){
  
  lim<-max(abs(x), na.rm=T)
  breaks<-seq(from=-1*lim, to=lim, length.out=n+1)
  class<-rep(NA, length(x))
  for(ii in 1:n){
    class[x>=breaks[ii] & x<breaks[ii+1]]<-ii
  }
  return(list(class=class, breaks=breaks))
  
}


legtext<-function(breaks, digits=2){
  out<-NULL
  for(ii in 1:(length(breaks)-1)){
    out<-c(out, paste0(round(breaks[ii],digits)," to ", round(breaks[ii+1],digits)))
  }
  return(out)
}


# pdf("sync_taildiff_25km.pdf")
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(coravg.diff.25km,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="mean <=25 km synchrony lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(coravg.diff.25km,nclasses)$breaks, digits=3),
#        col=pal)
# dev.off()


# pdf("sync_taildiff_25km_3panel.pdf", width=6.5, height=3)
# 
# par(mfrow=c(1,3),mar=c(4.1,4.1,2,1))
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(coravg.lb.25km,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(coravg.lb.25km,nclasses)$breaks, digits=3)
#        , col=pal, cex=0.8)
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(coravg.ub.25km,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(coravg.ub.25km,nclasses)$breaks, digits=3), 
#        col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(coravg.diff.25km,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(coravg.diff.25km,nclasses)$breaks, digits=3)
#        , col=pal, cex=0.8)
# 
# dev.off()


#Tail dependence in environmental relationships

kelpXwaves.all<-rep(NA, nrow(kelp))
kelpXwaves.lb<-rep(NA, nrow(kelp))
kelpXwaves.ub<-rep(NA, nrow(kelp))
kelpXno3.all<-rep(NA, nrow(kelp))
kelpXno3.lb<-rep(NA, nrow(kelp))
kelpXno3.ub<-rep(NA, nrow(kelp))
kelpXnpgo.lb<-rep(NA, nrow(kelp))
kelpXnpgo.ub<-rep(NA, nrow(kelp))
kelpXnpgo.all<-rep(NA, nrow(kelp))

for(ii in 1:nrow(kelp)){
  
  kelpXwaves.lb[ii]<-partialSpearman(kelp[ii,],-1*waves[ii,],c(0,0.5))
  kelpXwaves.ub[ii]<-partialSpearman(kelp[ii,],-1*waves[ii,],c(0.5,1))
  kelpXwaves.all[ii]<-cor(kelp[ii,],-1*waves[ii,], method="spearman")
  kelpXno3.lb[ii]<-partialSpearman(kelp[ii,],no3[ii,],c(0,0.5))
  kelpXno3.ub[ii]<-partialSpearman(kelp[ii,],no3[ii,],c(0.5,1))
  kelpXno3.all[ii]<-cor(kelp[ii,],no3[ii,], method="spearman")
  kelpXnpgo.lb[ii]<-partialSpearman(kelp[ii,],clims$NPGO,c(0,0.5))
  kelpXnpgo.ub[ii]<-partialSpearman(kelp[ii,],clims$NPGO,c(0.5,1))
  kelpXnpgo.all[ii]<-cor(kelp[ii,],clims$NPGO, method="spearman")
  
}


kelpXwaves.filt<-ifelse(kelpXwaves.all>=0,1,NA)
kelpXno3.filt<-ifelse(kelpXno3.all>=0,1,NA)
kelpXnpgo.filt<-ifelse(kelpXnpgo.all>=0,1,NA)

kelpXwaves.diff<-kelpXwaves.ub-kelpXwaves.lb*kelpXwaves.filt
kelpXno3.diff<-kelpXno3.ub-kelpXno3.lb*kelpXno3.filt
kelpXnpgo.diff<-kelpXnpgo.ub-kelpXnpgo.lb*kelpXnpgo.filt

# pdf("hist_taildiff_kelpXwaves.pdf")
# hist(kelpXwaves.diff, xlab="lower tail - upper tail")
# abline(v=median(kelpXwaves.diff, na.rm=T), col="red", lty=2)
# dev.off()

# pdf("hist_taildiff_kelpXno3.pdf")
# hist(kelpXno3.diff, xlab="lower tail - upper tail")
# abline(v=median(kelpXno3.diff,na.rm=T), col="red", lty=2)
# dev.off()
# 
# pdf("hist_taildiff_kelpXnpgo.pdf")
# hist(kelpXnpgo.diff, xlab="lower tail - upper tail")
# abline(v=median(kelpXnpgo.diff,na.rm=T), col="red", lty=2)
# dev.off()

cor(kelpXwaves.diff, kelpXno3.diff, use="pairwise.complete.obs") #these are slightly positively correlated

cor(kelpXwaves.diff, kelpXnpgo.diff, use="pairwise.complete.obs") #these are moderately negatively correlated

cor(kelpXno3.diff, kelpXnpgo.diff, use="pairwise.complete.obs") #these are moderately positively correlated

# pdf("map_taildiff_kelpXwaves.pdf")
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXwaves.diff,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="Waves lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXwaves.diff,nclasses)$breaks),
#        col=pal)
# dev.off()


# pdf("map_taildiff_kelpXwaves_3panel.pdf", width=6.5, height=3)
# 
# par(mfrow=c(1,3),mar=c(4.1,4.1,2,1))
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXwaves.lb*kelpXwaves.filt,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXwaves.lb*kelpXwaves.filt,nclasses)$breaks), col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXwaves.ub*kelpXwaves.filt,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXwaves.ub*kelpXwaves.filt,nclasses)$breaks), col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXwaves.diff,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXwaves.diff*kelpXwaves.filt,nclasses)$breaks), col=pal, cex=0.8)
# dev.off()


# pdf("map_taildiff_kelpXno3.pdf")
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXno3.diff,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="NO3 lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXno3.diff,nclasses)$breaks),
#        col=pal)
# dev.off()


# pdf("map_taildiff_kelpXno3_3panel.pdf", width=6.5, height=3)
# 
# par(mfrow=c(1,3),mar=c(4.1,4.1,2,1))
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXno3.lb*kelpXno3.filt,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXno3.lb*kelpXno3.filt,nclasses)$breaks), col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXno3.ub*kelpXno3.filt,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXno3.ub*kelpXno3.filt,nclasses)$breaks), col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXno3.diff,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXno3.diff,nclasses)$breaks), col=pal, cex=0.8)
# dev.off()



# pdf("map_taildiff_kelpXnpgo.pdf")
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXnpgo.diff,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="NPGO lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXnpgo.diff,nclasses)$breaks),
#        col=pal)
# dev.off()



# pdf("map_taildiff_kelpXnpgo_3panel.pdf", width=6.5, height=3)
# 
# par(mfrow=c(1,3),mar=c(4.1,4.1,2,1))
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXnpgo.lb*kelpXnpgo.filt,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXnpgo.lb*kelpXnpgo.filt,nclasses)$breaks), col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXnpgo.ub*kelpXnpgo.filt,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXnpgo.ub*kelpXnpgo.filt,nclasses)$breaks), col=pal, cex=0.8)
# 
# plot(locs$Lon, locs$Lat, pch=19, col=pal[class_eq_int_sym(kelpXnpgo.diff,nclasses)$class], cex=0.5,
#      xlab="Longitude", ylab="Latitude", main="lower tail - upper tail")
# legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXno3.diff,nclasses)$breaks), col=pal, cex=0.8)
# dev.off()


#make the environmental covariates averages within 50km!

avg_by_dist2<-function(x, dmat, thresh){
  
  nbhd<-ifelse(dmat <= thresh, 1, 0)
  nbhd.cor<-x * nbhd
  return(apply(nbhd.cor, 1, mean, na.rm=T))
  
}




moddat<-data.frame(kelp = coravg.diff.25km,
                   no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 25),
                   waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 25),
                   npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 25),
                   lat = locs$Lat,
                   lon = locs$Lon)

moddat<-moddat[complete.cases(moddat),]

mod<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat, nugget=T), data=moddat)
summary(mod)

rescor.mod<-spline.correlog(x=moddat$lon, y=moddat$lat, z=resid(mod), latlon=TRUE)
qq.mod<-qqnorm(mod)

png("./manuscript/diag1.png", res=300, units="in", width=6.5, height=3)

par(mar=c(3.1,3.1,1.1,1.1), mfrow=c(1,2), mgp=c(2,0.8,0))

plot(qq.mod$panel.args[[1]]$x, qq.mod$panel.args[[1]]$y,
     xlab="Standardized residuals",
     ylab="Quantiles of standard normal", pch=16, cex=0.7)
abline(a=0,b=1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")
plot(rescor.mod, xlab="Distance (km)")
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"b)")

dev.off()


#calculate spatial similarity following Bjornstad et al. 


resdist <- dist(resid(mod))

moddist <- gcdist(moddat$lon, moddat$lon)

plot(moddist[lower.tri(moddist)], resdist, pch=16, cex=0.7)


## Test sensitivity to distance threshold


## distance =  10 km
moddat3<-data.frame(kelp = avg_by_dist(cormat.diff, geog.dist, 10),
                    no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 10),
                    waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 10),
                    npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 10),
                    lat = locs$Lat,
                    lon = locs$Lon)

moddat3<-moddat3[complete.cases(moddat3),]

mod3<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat, nugget=T), data=moddat3)
summary(mod3)

## distance = 50 km
moddat4<-data.frame(kelp = avg_by_dist(cormat.diff, geog.dist, 50),
                    no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 50),
                    waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 50),
                    npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 50),
                    lat = locs$Lat,
                    lon = locs$Lon)

moddat4<-moddat4[complete.cases(moddat4),]

mod4<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat, nugget=T), data=moddat4)
summary(mod4)

## distance = 100 km
moddat5<-data.frame(kelp = avg_by_dist(cormat.diff, geog.dist, 100),
                    no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 100),
                    waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 100),
                    npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 100),
                    lat = locs$Lat,
                    lon = locs$Lon)

moddat5<-moddat5[complete.cases(moddat5),]

mod5<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat, nugget=T), data=moddat5)
summary(mod5)

## distance = 200 km
moddat6<-data.frame(kelp = avg_by_dist(cormat.diff, geog.dist, 200),
                    no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 200),
                    waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 200),
                    npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 200),
                    lat = locs$Lat,
                    lon = locs$Lon)

moddat6<-moddat6[complete.cases(moddat6),]

mod6<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat, nugget=T), data=moddat6)
summary(mod6)

#
# pdf("scatterplot_by_association.pdf")
# par(xpd=T)
# plot(NA, NA, xlim=range(-1*waves.raw), ylim=range(kelp.raw), xlab="calmness", ylab="kelp biomass")
# legend('top', pch=16, col=c("blue","red"),legend=c("lower stronger","upper stronger"), ncol=2, inset=-0.08)
# 
# for(ii in 1:nrow(waves)){
#   
#   if(is.na(kelpXwaves.diff[ii])){next}
#   
#   else if(kelpXwaves.diff[ii] < 0){
#     points(-1*waves.raw[ii,], kelp.raw[ii,], pch=16, col="blue")
#   }
#   else if(kelpXwaves.diff[ii] > 0){
#     points(-1*waves.raw[ii,], kelp.raw[ii,], pch=16, col="red")
#   }
#   
# }
# 
# dev.off()


# pdf("example_strong_upper.pdf")
# par(mfrow=c(2,1), mar=c(4.1,4.1,2.1,4.1))
# plot(years, kelp[which.min(kelpXwaves.diff),], type="b", col="blue", ylab="Kelp biomass (residuals)",
#      main="strong upper tail association")
# par(new=T)
# plot(years, -1*waves[which.min(kelpXwaves.diff),], type="b", col="red", yaxt="n", ylab="")
# plot(-1*waves[which.min(kelpXwaves.diff),], kelp[which.min(kelpXwaves.diff),],
#      xlab="calmness (residuals)", ylab="kelp biomass (residuals)")
# dev.off()
# 
# pdf("example_strong_lower.pdf")
# par(mfrow=c(2,1), mar=c(4.1,4.1,2.1,4.1))
# plot(years, kelp[which.max(kelpXwaves.diff),], type="b", col="blue", ylab="Kelp biomass (residuals)",
#      main="strong lower tail association")
# par(new=T)
# plot(years, -1*waves[which.max(kelpXwaves.diff),], type="b", col="red", yaxt="n", ylab="")
# plot(-1*waves[which.max(kelpXwaves.diff),], kelp[which.max(kelpXwaves.diff),],
#      xlab="calmness (residuals)", ylab="kelp biomass (residuals)")
# dev.off()


cor.test(rowMeans(kelp.raw), kelpXwaves.diff, method="spearman")

cor.test(-1*rowMeans(waves.raw), kelpXwaves.diff)

moddf<-data.frame(tail = kelpXwaves.diff,
                  calm = -1*rowMeans(waves.raw),
                  lat = locs$Lat,
                  lon = locs$Lon)
moddf<-moddf[complete.cases(moddf),]

mod_tail_v_calm <- gls(tail ~ calm, correlation = corExp(form = ~lat + lon), data=moddf)
summary(mod_tail_v_calm)
qq.mod<-qqnorm(mod_tail_v_calm)

rescor.mod<-spline.correlog(moddf$lon, moddf$lat, resid(mod_tail_v_calm), latlon=TRUE)

png("./manuscript/diag2.png", res=300, units="in", width=6.5, height=3)

par(mar=c(3.1,3.1,1.1,1.1), mfrow=c(1,2), mgp=c(2,0.8,0))

plot(qq.mod$panel.args[[1]]$x, qq.mod$panel.args[[1]]$y,
     xlab="Standardized residuals",
     ylab="Quantiles of standard normal", pch=16, cex=0.7)
abline(a=0,b=1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")
plot(rescor.mod, xlab="Distance (km)")
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"b)")

dev.off()


# pdf("calmness_vs_tail.pdf")
# plot(-1*rowMeans(waves.raw), kelpXwaves.diff, ylab="kelpXwaves lower - upper", xlab="Mean wave calmness")
# dev.off()



## Investigate skewness of time series


sample_skewness<-function(x){
  #from Joanes, D. N.; Gill, C. A. (1998). "Comparing measures of sample skewness and kurtosis". 
  # Journal of the Royal Statistical Society, Series D. 47 (1): 183–189. doi:10.1111/1467-9884.00122.
  #this form is shown to have low RMSE in small samples from non-normal distributions
  if(any(is.na(x))){
    stop("x contains one or more NA values")
  }
  
  n = length(x)
  xbar = mean(x)
  m3 = (1/n)*sum((x-xbar)^3)
  s3 = ((1/(n-1))*sum((x-xbar)^2))^(3/2)
  c = sqrt(n*(n-1))/(n-2)
  return(c*m3/s3)
  
}


## this is for individual kelp time series
kelp_skew<-apply(kelp.raw, 1, sample_skewness)
hist(kelp_skew) #distributions are primarily right-skewed (high outliers)


plot(coravg.diff.25km, kelp_skew) ## <-------------- this needs to be for an aggregate!
cor.test(coravg.diff.25km, kelp_skew, method="spearman")


moddat2<-data.frame(kelp_skew = kelp_skew,
                   no3 = kelpXno3.diff,
                   waves = kelpXwaves.diff,
                   npgo =kelpXnpgo.diff,
                   lat = locs$Lat,
                   lon = locs$Lon)

moddat2<-moddat2[complete.cases(moddat2),]

mod2<-gls(kelp_skew ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat, nugget=T), data=moddat2)
summary(mod2) #no statistical significant effects


#this is for spatially averaged kelp time series
dthresh = 25 #km
neighmat<-ifelse(geog.dist < dthresh, 1, 0)

kelp_spatavg<-matrix(NA, nrow(kelp), ncol(kelp))

for(ii in 1:nrow(kelp.raw)){
  for(jj in 1:ncol(kelp.raw)){
    
    kelp_spatavg[ii,jj] <- mean(kelp.raw[,jj]*neighmat[ii,])
    
  }
}

kelp_spatavg_skew<-apply(kelp_spatavg, 1, sample_skewness)

hist(kelp_spatavg_skew) #distributions are primarily right-skewed (high outliers)

plot(coravg.diff.25km, kelp_spatavg_skew, xlim=c(-0.03,0.03))
abline(v=0, col="red", lty=2)
abline(h=0, col="blue", lty=2)
mtext("lower tail association", at=-0.02)
mtext("upper tail association", at=0.02)
cor.test(coravg.diff.25km, kelp_spatavg_skew, method="spearman") 
# positive correlation but funky dist'n so using Spearman; rho = 0.27, p << 0.0001

cor.test(coravg.diff.25km, rowMeans(kelp_spatavg))


moddat3<-data.frame(kelp = kelp_spatavg_skew,
                   no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 25),
                   waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 25),
                   npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 25),
                   lat = locs$Lat,
                   lon = locs$Lon)

moddat3<-moddat3[complete.cases(moddat3),]

mod3<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(value=0.01, form = ~ lon + lat, nugget=T), data=moddat3)
summary(mod3) #no significant effects; model does not fit without setting value of corExp, but results are
# consistent across a range of effects


## Make manuscript figures ------------------------------------------------------------------------

#Fig 1: distance decay

# png("./manuscript/fig2_distdecay.png", width=6.5, height=3, units="in", res=300)
# 
# par(mfrow=c(1,3), mgp=c(2.3,0.7,0), mar=c(2.1,2.1,1.6,1.1), cex.lab=1.3, oma=c(2,2,0,0))
# 
# plot(sncf.lb$emp.spline$x, sncf.lb$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
#      xlab="", ylab="", main="All coastline segments")
# lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[1,], col="red", lty=1)
# lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[3,], col="red", lty=1)
# lines(sncf.ub$emp.spline$x, sncf.ub$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[1,], col="blue", lty=1)
# lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[3,], col="blue", lty=1)
# text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")
# legend("topright",lty=1,col=c("red","blue"),legend=c("lower tail","upper tail"),bty="n")
# 
# plot(sncf.lb.c1$emp.spline$x, sncf.lb.c1$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
#      xlab="", ylab="", main="Central California")
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[1,], col="red", lty=1)
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[3,], col="red", lty=1)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[1,], col="blue", lty=1)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[3,], col="blue", lty=1)
# text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"b)")
# #legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
# 
# plot(sncf.lb.c2$emp.spline$x, sncf.lb.c2$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
#      xlab="", ylab="", main="Southern California")
# lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[1,], col="red", lty=1)
# lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[3,], col="red", lty=1)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[1,], col="blue", lty=1)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[3,], col="blue", lty=1)
# text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"c)")
# #legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
# 
# mtext("Distance (km)", side=1, outer=T, line=0.3)
# mtext("Partial Spearman correlation", side=2, outer=T, line=0.2)
# 
# dev.off()


## Fig 2: Matrices

pal<-colorRampPalette(colors=c("red","white","blue"))

axis.at=c(32,75,211,243,315,359)
axis.labels=c("Monterey","Point Sur","Morro Bay","Pt. Conception","Santa Monica","San Diego")
axis.labels2<-c("MO","PS","MB","PC","SM","SD")


## NEW FIGURE 3 COMBINES DISTANCE DECAY AND MATRICES

laymat <- matrix(1:9, nrow=3, ncol=3)
wr <- 0.07

png("./manuscript/fig3_combSynchrony.png", width=6.25, height=8.5, units="in", res=300)

layout(laymat, widths=c((1-wr)/2, (1-wr)/2, wr))

par(mgp=c(2.3,0.7,0), mar=c(4.1,4.1,1.6,1.1), cex.lab=1.2, oma=c(6,0,0,0))

plot(sncf.lb$emp.spline$x, sncf.lb$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
     xlab="Distance", ylab="Partial Spearman correlation", main="")
lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[1,], col="red", lty=1)
lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[3,], col="red", lty=1)
lines(sncf.ub$emp.spline$x, sncf.ub$emp.spline$y, col="blue", lwd=2)
lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[1,], col="blue", lty=1)
lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[3,], col="blue", lty=1)
mtext("a)", cex=0.75, at=par("usr")[1], line=0.2)
mtext("All coastline segments", cex=0.7, line=0.2)
#text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")
legend("topright",lty=1,col=c("red","blue"),legend=c("lower tail","upper tail"),bty="n")

plot(sncf.lb.c1$emp.spline$x, sncf.lb.c1$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
     xlab="Distance", ylab="Partial Spearman correlation", main="")
lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[1,], col="red", lty=1)
lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[3,], col="red", lty=1)
lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$emp.spline$y, col="blue", lwd=2)
lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[1,], col="blue", lty=1)
lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[3,], col="blue", lty=1)
#text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"b)")
mtext("b)", cex=0.75, at=par("usr")[1], line=0.2)
mtext("Central California", cex=0.7, line=0.2)
#legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))

plot(sncf.lb.c2$emp.spline$x, sncf.lb.c2$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
     xlab="Distance", ylab="Partial Spearman correlation", main="")
lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[1,], col="red", lty=1)
lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[3,], col="red", lty=1)
lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$emp.spline$y, col="blue", lwd=2)
lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[1,], col="blue", lty=1)
lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[3,], col="blue", lty=1)
#text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"c)")
mtext("c)", cex=0.75, at=par("usr")[1], line=0.2)
mtext("Southern California", cex=0.7, line=0.2)
#legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))


par(mgp=c(2,0.7,0), mar=c(1.5,2.1,1.6,1.1))

#image.plot(1:361,1:361,cormat.ub.filt, xaxt="n", yaxt="n", legend.width=0.75, legend.mar=2.5,
#           zlim=c(-0.6,0.6), col=pal(64))
image(1:361,1:361,cormat.ub.filt, xaxt="n", yaxt="n",
           zlim=c(-0.6,0.6), col=pal(64))
mtext("d)", cex=0.75, at=0, line=0.2)
mtext("Upper tail synchrony matrix", cex=0.7, line=0.2)
axis(side=1,at=axis.at, labels=axis.labels2)
axis(side=2,at=axis.at, labels=axis.labels2, las=2)
# image.plot(1:361,1:361,cormat.lb.filt, xaxt="n", yaxt="n", legend.width=0.75, legend.mar=2.5,
#            zlim=c(-0.6,0.6), col=pal(64))
image(1:361,1:361,cormat.lb.filt, xaxt="n", yaxt="n",
           zlim=c(-0.6,0.6), col=pal(64))
mtext("e)", cex=0.75, at=0, line=0.2)
mtext("Lower tail synchrony matrix", cex=0.7, line=0.2)
axis(side=1,at=axis.at, labels=axis.labels2)
axis(side=2,at=axis.at, labels=axis.labels2, las=2)
# image.plot(1:361,1:361,cormat.diff.nodiag, xaxt="n", yaxt="n", legend.width=0.75, legend.mar=2.5,
#            zlim=c(-0.6,0.6), col=pal(64))
image(1:361,1:361,cormat.diff.nodiag, xaxt="n", yaxt="n",
           zlim=c(-0.6,0.6), col=pal(64))
mtext("f)", cex=0.75, at=0, line=0.2)
mtext("Tail dependence strength (upper-lower)", cex=0.7, line=0.2)
axis(side=1,at=axis.at, labels=axis.labels, las=2)
axis(side=2,at=axis.at, labels=axis.labels2, las=2)


par(mar=c(1.5,1.5,1.6,1.1))
image(t(matrix(1:64)), xaxt="n", yaxt="n", col=pal(64))
axis(2, at=seq(0,1,length.out=5),labels=seq(-0.6,0.6,length.out=5))
image(t(matrix(1:64)), xaxt="n", yaxt="n", col=pal(64))
axis(2, at=seq(0,1,length.out=5),labels=seq(-0.6,0.6,length.out=5))
image(t(matrix(1:64)), xaxt="n", yaxt="n", col=pal(64))
axis(2, at=seq(0,1,length.out=5),labels=seq(-0.6,0.6,length.out=5))

dev.off()




# png("./manuscript/fig3_matrices.png", width=3.2, height=8.5, units="in", res=300)
# 
# par(mfrow=c(3,1), mgp=c(2,0.5,0), mar=c(1.1,1.1,1.5,1.1), oma=c(6,1.5,0,0))
# 
# image.plot(1:361,1:361,cormat.ub.filt, xaxt="n", yaxt="n", legend.width=0.75, legend.mar=2.5, 
#            zlim=c(-0.6,0.6), col=pal(64))
# mtext("a)", cex=0.75, at=0, line=0.2)
# mtext("Upper tail synchrony matrix", cex=0.7)
# axis(side=1,at=axis.at, labels=axis.labels2)
# axis(side=2,at=axis.at, labels=axis.labels2, las=2)
# image.plot(1:361,1:361,cormat.lb.filt, xaxt="n", yaxt="n", legend.width=0.75, legend.mar=2.5, 
#            zlim=c(-0.6,0.6), col=pal(64))
# mtext("b)", cex=0.75, at=0, line=0.2)
# mtext("Lower tail synchrony matrix", cex=0.7)
# axis(side=1,at=axis.at, labels=axis.labels2)
# axis(side=2,at=axis.at, labels=axis.labels2, las=2)
# image.plot(1:361,1:361,cormat.diff.nodiag, xaxt="n", yaxt="n", legend.width=0.75, legend.mar=2.5, 
#            zlim=c(-0.6,0.6), col=pal(64))
# mtext("c)", cex=0.75, at=0, line=0.2)
# mtext("Tail dependence strength (upper-lower)", cex=0.7)
# axis(side=1,at=axis.at, labels=axis.labels, las=2)
# axis(side=2,at=axis.at, labels=axis.labels2, las=2)
# 
# dev.off()


## Figure 4: Kelp synchrony map

#states<-readOGR("/Users/jonathanwalter/Documents/Research/DATA/Basemaps/tl_us_states_2000/tl_2009_us_state00.shp")
#cali<-states[states$NAME00=="California",]
# plot(cali)
# proj4string(cali)

for (package in c('dplyr', 'tidyr', 'ggplot2', 'PBSmapping', 'animation', 'lme4', 'car', 'DHARMa', 'ggspatial',
                  'sjmisc', 'geosphere', 'raster', 'sp', 'rgdal', #'SDMTools', 
                  'ncf', 'maptools', 'multcomp',
                  'gridExtra', 'ROCR', 'pROC', 'raster', 'rasterVis', 'reshape2', 'RColorBrewer', 'rgeos')) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

# Import coastline data from NOAA
coast.limits <- list(x = c(-125 + 360, -114 + 360), y = c(30, 39))

# Coastline polygons
coast.polys <- importGSHHS(gshhsDB="~/Box Sync/Coastline vectors from GSHHG/gshhg-bin-2.3.7/gshhs_f.b",
                           xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
coast.polys$X <- coast.polys$X - 360  # Get longitude back into units of degrees east
coast.polys.sp <- maptools::PolySet2SpatialPolygons(coast.polys, close_polys=FALSE) # Convert to spatial polygons
polys <- fortify(coast.polys)

# # Coastline borders
# coast.border <- importGSHHS(gshhsDB="~/Box Sync/Coastline vectors from GSHHG/gshhg-bin-2.3.7/wdb_borders_f.b",
#                             xlim=coast.limits$x, ylim=coast.limits$y, maxLevel=1, n=0)
# coast.border$X <- coast.border$X - 360  # Get longitude back into units of degrees east
# coast.border.sp <- maptools::PolySet2SpatialPolygons(coast.border, close_polys=FALSE) # Convert to spatial polygons
# border <- fortify(coast.border)

#plot(coast.border.sp)

# thin <- c(seq(1, 298, by=5),seq(299, nrow(locs), by=2))
thin <- c(seq(1,211,by=5),seq(211,243,by=2),seq(243,299,by=5),seq(299,nrow(locs),by=2))

places <- data.frame(x = c(-121.8, -121.7, -120.8, -120.4, -118.49, -117.15),
                     y = c(36.6, 36.3, 35.37, 34.55, 34.02, 32.72),
                     name = c("Monterey","Point Sur","Morro Bay", "Point Conception", "Santa Monica", "San Diego"))


nclasses=8
pal<-brewer.pal(nclasses,"RdYlBu")
tpal<-c("#D7302780", "#F46D4380", "#FDAE6180", "#FEE09080", "#E0F3F880", "#ABD9E980", "#74ADD180", "#4575B480")


png("./manuscript/fig4_kelpsynch_maps.png", width=6.5, height=3, res=300, units="in")

par(mfrow=c(1,3), mar=c(0.5,0.5,2,0.5))

plot(coast.polys.sp, xlim=range(locs$Lon), ylim=range(locs$Lat))
points(locs$Lon[thin], locs$Lat[thin], pch=19, 
       col=pal[class_eq_int_sym(coravg.ub.25km,nclasses)$class][thin], cex=0.6)
     #xlab="Longitude", ylab="Latitude", main="lower tail")
mtext("a)",at=par("usr")[1]+0.05*diff(par("usr")[1:2]), side=3, cex=0.75, line=0.5)
mtext("Upper tail synchrony", cex=0.75, line=0.5)
#rect(bbox(cali)[1],bbox(cali)[2],bbox(cali)[3],bbox(cali)[4])
legend("topright", pch=19, 
       legend=legtext(class_eq_int_sym(coravg.lb.25km,nclasses)$breaks, digits=3)[5:8]
       , col=pal[5:8], cex=0.8)
points(places$x, places$y, pch=23, col="grey", bg="black", lwd=0.5, cex=0.8)
text(places$x, places$y, places$name, pos=c(rep(4,5),2), cex=0.65)

plot(coast.polys.sp, xlim=range(locs$Lon), ylim=range(locs$Lat))
points(locs$Lon[thin], locs$Lat[thin], pch=19, 
       col=pal[class_eq_int_sym(coravg.lb.25km,nclasses)$class][thin], cex=0.6,
     xlab="Longitude", ylab="Latitude", main="upper tail")
mtext("Lower tail synchrony", cex=0.75, line=0.5)
mtext("b)",at=par("usr")[1]+0.05*diff(par("usr")[1:2]), side=3, cex=0.75, line=0.5)
legend("topright", pch=19, legend=legtext(class_eq_int_sym(coravg.ub.25km,nclasses)$breaks, digits=3)[5:8], 
       col=pal[5:8], cex=0.8)
points(places$x, places$y, pch=23, col="grey", bg="black", lwd=0.5, cex=0.8)
text(places$x, places$y, places$name, pos=c(rep(4,5),2), cex=0.65)

plot(coast.polys.sp, xlim=range(locs$Lon), ylim=range(locs$Lat))
points(locs$Lon[thin], locs$Lat[thin], pch=19, 
       col=pal[class_eq_int_sym(coravg.diff.25km,nclasses)$class][thin], cex=0.6,
     xlab="Longitude", ylab="Latitude", main="lower tail - upper tail")
mtext("c)",at=par("usr")[1]+0.05*diff(par("usr")[1:2]), side=3, cex=0.75, line=0.5)
mtext("Tail dependence strength", cex=0.75, line=0.5)
legend("topright", pch=19, legend=legtext(class_eq_int_sym(coravg.diff.25km,nclasses)$breaks, digits=3)
       , col=pal, cex=0.8)
segments(-121.8, 36.7, -120.85, 37.735, col="grey") #upper inset
segments(-121.45, 35.92, -120.85, 35.8, col="grey")
segments(-120.649, 34.65, -121.995, 33.38, col="grey") #lower inset
segments(-119.657, 34.355, -119.195, 33.38, col="grey")
points(places$x, places$y, pch=23, col="grey", bg="black", lwd=0.5, cex=0.8)
text(places$x, places$y, c("MO","PS",places$name[3:6]), pos=c(rep(4,5),2), cex=0.65)

par(new=T, fig=c(0.77, 0.85, 0.6, 0.9), mar=c(0.05,0.05,0.05,0.05)) #upper inset
plot(coast.polys.sp, ylim=c(36.1,36.5), xlim=c(-122,-121.5))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(coravg.diff.25km,nclasses)$class], cex=0.5)
rect(-122.005, 35.925, -121.485, 36.676, col=NA, border="grey")

par(new=T, fig=c(0.7, 0.87, 0.065, 0.24), mar=c(0.05,0.05,0.05,0.05)) #lower inset
plot(coast.polys.sp, ylim=c(34.3,34.6), xlim=c(-120.7,-119.7))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(coravg.diff.25km,nclasses)$class], cex=0.5)
rect(-120.7, 34.3, -119.663, 34.655, col=NA, border="grey")


dev.off()


## Figure 4: Driver tail association maps
png("./manuscript/fig5_drivertail_maps.png", width=6.5, height=3, res=300, units="in")

par(mfrow=c(1,3),mar=c(0.5,0.5,2,0.5))

plot(coast.polys.sp, xlim=range(locs$Lon), ylim=range(locs$Lat))
points(locs$Lon[thin], locs$Lat[thin], pch=19, 
       col=pal[class_eq_int_sym(kelpXwaves.diff,nclasses)$class][thin], cex=0.6,
       xlab="Longitude", ylab="Latitude", main="lower tail")
mtext("a)",at=par("usr")[1]+0.05*diff(par("usr")[1:2]), side=3, cex=0.75, line=0.5)
mtext("Wave calmness", cex=0.75, line=0.5)
#rect(bbox(cali)[1],bbox(cali)[2],bbox(cali)[3],bbox(cali)[4])
legend("topright", pch=19, 
       legend=legtext(class_eq_int_sym(kelpXwaves.diff,nclasses)$breaks, digits=3)
       , col=pal, cex=0.8)
segments(-121.8, 36.7, -120.85, 37.735, col="grey")
segments(-121.45, 35.92, -120.85, 35.8, col="grey")
segments(-120.649, 34.65, -121.995, 33.38, col="grey") #lower inset
segments(-119.657, 34.355, -119.195, 33.38, col="grey")
points(places$x, places$y, pch=23, col="grey", bg="black", lwd=0.5, cex=0.8)
text(places$x, places$y, c("MO","PS",places$name[3:6]), pos=c(rep(4,5),2), cex=0.65)

plot(coast.polys.sp, xlim=range(locs$Lon), ylim=range(locs$Lat))
points(locs$Lon[thin], locs$Lat[thin], pch=19, 
       col=pal[class_eq_int_sym(kelpXno3.diff,nclasses)$class][thin], cex=0.6,
       xlab="Longitude", ylab="Latitude", main="upper tail")
mtext("b)",at=par("usr")[1]+0.05*diff(par("usr")[1:2]), side=3, cex=0.75, line=0.5)
mtext("Nitrate concentration", cex=0.75, line=0.5)
legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXno3.diff,nclasses)$breaks, digits=3), 
       col=pal, cex=0.8)
segments(-121.8, 36.7, -120.85, 37.735, col="grey")
segments(-121.45, 35.92, -120.85, 35.8, col="grey")
segments(-120.649, 34.65, -121.995, 33.38, col="grey") #lower inset
segments(-119.657, 34.355, -119.195, 33.38, col="grey")
points(places$x, places$y, pch=23, col="grey", bg="black", lwd=0.5, cex=0.8)
text(places$x, places$y, c("MO","PS",places$name[3:6]), pos=c(rep(4,5),2), cex=0.65)

plot(coast.polys.sp, xlim=range(locs$Lon), ylim=range(locs$Lat))
points(locs$Lon[thin], locs$Lat[thin], pch=19, 
       col=pal[class_eq_int_sym(kelpXnpgo.diff,nclasses)$class][thin], cex=0.6,
       xlab="Longitude", ylab="Latitude", main="lower tail - upper tail")
mtext("c)",at=par("usr")[1]+0.05*diff(par("usr")[1:2]), side=3, cex=0.75, line=0.5)
mtext("NPGO", cex=0.75, line=0.5)
legend("topright", pch=19, legend=legtext(class_eq_int_sym(kelpXnpgo.diff,nclasses)$breaks, digits=3)
       , col=pal, cex=0.8)
segments(-121.8, 36.7, -120.85, 37.735, col="grey")
segments(-121.45, 35.92, -120.85, 35.8, col="grey")
segments(-120.649, 34.65, -121.995, 33.38, col="grey") #lower inset
segments(-119.657, 34.355, -119.195, 33.38, col="grey")
points(places$x, places$y, pch=23, col="grey", bg="black", lwd=0.5, cex=0.8)
text(places$x, places$y, c("MO","PS",places$name[3:6]), pos=c(rep(4,5),2), cex=0.65)

par(new=T, fig=c(0.104, 0.184, 0.6, 0.9), mar=c(0.05,0.05,0.05,0.05))
plot(coast.polys.sp, ylim=c(36.1,36.5), xlim=c(-122,-121.5))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(kelpXwaves.diff,nclasses)$class], cex=0.5)
rect(-122.005, 35.925, -121.485, 36.676, col=NA, border="grey")

par(new=T, fig=c(0.437, 0.517, 0.6, 0.9), mar=c(0.05,0.05,0.05,0.05))
plot(coast.polys.sp, ylim=c(36.1,36.5), xlim=c(-122,-121.5))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(kelpXno3.diff,nclasses)$class], cex=0.5)
rect(-122.005, 35.925, -121.485, 36.676, col=NA, border="grey")

par(new=T, fig=c(0.77, 0.85, 0.6, 0.9), mar=c(0.05,0.05,0.05,0.05))
plot(coast.polys.sp, ylim=c(36.1,36.5), xlim=c(-122,-121.5))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(kelpXnpgo.diff,nclasses)$class], cex=0.5)
rect(-122.005, 35.925, -121.485, 36.676, col=NA, border="grey")

par(new=T, fig=c(0.033, 0.203, 0.065, 0.24), mar=c(0.05,0.05,0.05,0.05)) #lower inset
plot(coast.polys.sp, ylim=c(34.3,34.6), xlim=c(-120.7,-119.7))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(kelpXwaves.diff,nclasses)$class], cex=0.5)
rect(-120.7, 34.3, -119.663, 34.655, col=NA, border="grey")

par(new=T, fig=c(0.366, 0.536, 0.065, 0.24), mar=c(0.05,0.05,0.05,0.05)) #lower inset
plot(coast.polys.sp, ylim=c(34.3,34.6), xlim=c(-120.7,-119.7))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(kelpXno3.diff,nclasses)$class], cex=0.5)
rect(-120.7, 34.3, -119.663, 34.655, col=NA, border="grey")

par(new=T, fig=c(0.7, 0.87, 0.065, 0.24), mar=c(0.05,0.05,0.05,0.05)) #lower inset
plot(coast.polys.sp, ylim=c(34.3,34.6), xlim=c(-120.7,-119.7))
points(locs$Lon, locs$Lat, pch=19,
       col=pal[class_eq_int_sym(kelpXnpgo.diff,nclasses)$class], cex=0.5)
rect(-120.7, 34.3, -119.663, 34.655, col=NA, border="grey")

dev.off()



## Figure 5: Calmness versus tail association


#schematic 

x <- seq(-6,6,by=0.1)
S <- 1/(1 + exp(-1.5*x))

# source("./code/extremeTailDep.R")
# np<-50
# xtc_r<-retd(n=np,d=2,rl=1,mn=0,sdev=1) # copula with extreme right tail
# 
# 
# # flip xtc_r to get xtc_l
# #xtc_l<-1-xtc_r
# # make marginals uniform 
# xtc_r<-pnorm(xtc_r)
# xtc_l<-1-xtc_r

## find strongly tail dependent relationships

find.strong <- cormat.diff.nodiag * ifelse(cormat.all > 0.5, 1, NA)
head(sort(find.strong), 10)
tail(sort(find.strong), 10)

which(find.strong < -0.3, arr.ind = TRUE)
which(find.strong > 0.3, arr.ind = TRUE)

which(find.strong==max(find.strong, na.rm=T), arr.ind = TRUE)


strong.waves <- kelpXwaves.diff * ifelse(kelpXwaves.all > 0.4, 1, NA)

head(sort(strong.waves),10)
tail(sort(strong.waves),10)

which(strong.waves < -0.15)
which(strong.waves > 0.15)



## Make figures for manuscript

png("./manuscript/fig6_calmnTailAssoc.png",units="in",res=300,height=6.5,width=6.5)

layout(matrix(c(1,2,3,3,4,4), nrow=3, ncol=2, byrow=TRUE), heights = c(0.5,0.3,0.3))

par(mar=c(3.1,3.1,3.1,1.1), mgp=c(2,0.8,0), cex.axis=0.9)

plot(-1*rowMeans(waves.raw), kelpXwaves.diff, ylab="Tail dependence (upper - lower)", 
     xlab="Mean wave calmness", pch=19)
abline(h=0, col="grey", lty=2)
abline(v=median(-1*rowMeans(waves.raw)), col="grey", lty=2)
abline(mod_tail_v_calm, col="green")
axis(3,at=seq(-6,-1,1),labels=rev(1:6))
mtext("Mean wave height", line=2, cex=0.65)
text(-4.5,-0.25,expression(paste(beta, "= -0.025, p = 0.001")), cex=1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")

plot(x, S, type="l", xaxt="n", yaxt="n", xlab="Wave calmness", ylab="Proportional kelp abundance",
     ylim=c(-0.22,1.22))
axis(1, at=seq(min(x),max(x),length.out=5), labels=c("Low","","","","High"))
axis(2, at=seq(0,1,length.out=5),labels=TRUE)
axis(3, at=seq(min(x),max(x),length.out=5), labels=rev(c("Low","","","","High")))
arrows(x0=-5, y0=-0.035, x1=-1, y1=-0.035, col="blue", length=0.03, angle=90, code=3, lwd=1.3)
text(-3,-.16,"Wave exposed site:\neffects in upper tail",col="blue")
text(3,1.16,"Calm site:\neffects in lower tail",col="red")
mtext("Wave height", line=2, cex=0.65)
arrows(x0=1, y0=1.035, x1=5, y1=1.035, col="red", length=0.03, angle=90, code=3, lwd=1.3)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"b)")

par(mar=c(3.1,3.1,1.25,1.1), mgp=c(2,0.8,0), cex.axis=0.9)

plot(NA, NA, xlim=c(1987,2019), ylim=c(-2,3), xlab="Time", ylab="Standardized kelp biomass")
lines(years, scale(kelp[48,]), col="blue")
lines(years, scale(kelp[54,]), lty=2, col="blue")
mtext("Upper tail association", cex=0.75, line=0.1)
text(par("usr")[1]+0.025*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"c)")
segments(x0=2005, x1=2009, y0=1.75, y1=2.3, col="darkgrey")
segments(x0=2013, x1=2009, y0=2.2, y1=2.3, col="darkgrey")
text(x=2009, y=2.3, "Synchronous booms", col="darkgrey", pos=3)
segments(x0=1998, x1=2002, y0=-1.9, y1=-1.6, col="darkgrey")
text(x=2002, y=-1.6, "Asynchronous crash", pos=4, col="darkgrey")


plot(NA, NA, xlim=c(1987,2019), ylim=c(-2,3), xlab="Time", ylab="Standardized kelp biomass")
lines(years, scale(kelp[30,]), col="red")
lines(years, scale(kelp[31,]), lty=2, col="red")
mtext("Lower tail association", cex=0.75, line=0.1)
text(par("usr")[1]+0.025*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"d)")
segments(x0=2000, x1=2005.5, y0=1.9, y1=2.3, col="darkgrey")
segments(x0=2018, x1=2012.5, y0=3, y1=2.6, col="darkgrey")
text(x=2009, y=2.45, "Asynchronous booms", col="darkgrey")
segments(x0=1992, x1=2002, y0=-1.3, y1=-1.6, col="darkgrey")
segments(x0=1995, x1=2002, y0=-1.6, y1=-1.6, col="darkgrey")
text(x=2002, y=-1.6, "Synchronous crashes", pos=4, col="darkgrey")

dev.off()



##

png("./manuscript/figSX_tailExamples.png",units="in",res=300,height=6.5,width=6.5)

par(mfrow=c(2,2), mar=c(3.5,3.5,1.6,1.1), mgp=c(2,0.7,0))

plot(rank(kelp[326,]), rank(kelp[330,]), xlab="Kelp biomass rank, site 326",
     ylab="Kelp biomass rank, site 330", pch=16, col="blue") #this looks OK for upper tail dependence
mtext("Upper tail dependence", cex=0.8, line=0.1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")

plot(rank(kelp[30,]), rank(kelp[28,]), xlab="Kelp biomass rank, site 30", 
     ylab="Kelp biomass rank, site 28", pch=16, col="red") #this looks pretty good for lower tail dependence
mtext("Lower tail dependence", cex=0.8, line=0.1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"b)")

plot(rank(-1*waves[44,]), rank(kelp[44,]), xlab="Wave calmness rank, site 44",
     ylab="Kelp biomass rank, site 44", pch=16, col="blue") #this looks OK for upper tail dependence
mtext("Upper tail dependence", cex=0.8, line=0.1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"c)")

plot(rank(-1*waves[245,]), rank(kelp[245,]), xlab="Wave calmess rank, site 245",
     ylab="Kelp biomass rank, site 245", pch=16, col="red") #this looks pretty good for lower tail dependence
mtext("Lower tail dependence", cex=0.8, line=0.1)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"d)")



dev.off()


