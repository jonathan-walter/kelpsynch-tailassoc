rm(list=ls())

library(wsyn)
library(ncf)
library(ecodist)
library(RColorBrewer)
library(nlme)
library(fields)

setwd("/Users/jonathanwalter/GitHub/kelp-synchrony/")

source("./TailAssociation/partialSpearman.R")

kelp<-as.matrix(read.csv("./CleanedData/Minimal/kelpbio.csv", header=F))
no3<-as.matrix(read.csv("./CleanedData/Minimal/NO3.csv", header=F))
waves<-as.matrix(read.csv("./CleanedData/Minimal/wavht.csv",header=F))
clims<-read.csv("./CleanedData/Minimal/clims.csv")
quarters<-read.csv("./CleanedData/Minimal/quarters.csv")
clusters<-read.csv("./CleanedData/Minimal/AllClustersNumeric.csv")

pdf("./TailAssociation/locs_by_cluster.pdf")
plot(clusters$Lon, clusters$Lat, col=clusters$kelpbio, pch=19, cex=0.2)
legend("topright",pch=19,col=1:2,legend=c("c1","c2"))
dev.off()

#make annual average time series
years<-unique(quarters$Year)
kelp.ann<-matrix(NA, nrow(kelp), length(years))
waves.ann<-matrix(NA, nrow(waves), length(years))
no3.ann<-matrix(NA, nrow(no3), length(years))
npgo.ann<-rep(NA, length(years))

for(ii in 1:nrow(kelp)){
  for(jj in 1:length(years)){
    kelp.ann[ii,jj]<-mean(kelp[ii,quarters$Year==years[jj]])
    waves.ann[ii,jj]<-mean(waves[ii,quarters$Year==years[jj]])
    no3.ann[ii,jj]<-mean(no3[ii,quarters$Year==years[jj]])
    if(ii==1){
      npgo.ann[jj]<-mean(clims$NPGO[quarters$Year==years[jj]])
    }
  }
}

#select persistent patches
persistent<-rowSums(kelp.ann==0)==0
kelp.ann<-kelp.ann[persistent,]
clusters<-clusters[persistent,]
waves.ann<-waves.ann[persistent,]
no3.ann<-no3.ann[persistent,]


lb<-c(0,0.5)
ub<-c(0.5,1)

cormat.all<-cor(t(kelp.ann), method="spearman")
filt.pos<-ifelse(cormat.all>=0,1,NA)

cormat.lb<-matrix(0,nrow(kelp.ann),nrow(kelp.ann))
cormat.ub<-matrix(0,nrow(kelp.ann),nrow(kelp.ann))
for(ii in 2:nrow(kelp.ann)){
  for(jj in 1:(ii-1)){
    
    cormat.lb[ii,jj]<-partialSpearman(kelp.ann[ii,],kelp.ann[jj,],lb)
    cormat.ub[ii,jj]<-partialSpearman(kelp.ann[ii,],kelp.ann[jj,],ub)
    
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
    shuffle<-sample(1:nrow(distmat), size=nrow(distmat), replace=FALSE)
    xres<-distmat[shuffle,shuffle][triang]
    yres<-zmat[shuffle,shuffle][triang]
    drop.NaNs<-!is.na(yres)
    xres<-xres[drop.NaNs]; yres<-yres[drop.NaNs]
    res.spline<-smooth.spline(xres,yres,df=dfs)
    resamp.splines[ii,]<-predict(res.spline, x=emp.spline$x)$y
  }
  out$resamp.splines<-resamp.splines
  out$spline.quantiles<-apply(resamp.splines,2,quantile,probs=quantiles)
  return(out)
}

geog.dist<-gcdist(clusters$Lon,clusters$Lat)

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

pdf("./TailAssociation/dist_decay_allsites.pdf")
plot(sncf.lb$emp.spline$x, sncf.lb$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
     xlab="Distance (km)", ylab="Partial Spearman correlation", main="All coastline segments")
# lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[1,], col="red", lty=3)
# lines(sncf.lb$emp.spline$x, sncf.lb$spline.quantiles[3,], col="red", lty=3)
lines(sncf.ub$emp.spline$x, sncf.ub$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[1,], col="blue", lty=3)
# lines(sncf.ub$emp.spline$x, sncf.ub$spline.quantiles[3,], col="blue", lty=3)
legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
dev.off()



#look at clusters separately -------------------------------------------------------------------------------------------

geog.dist.c1<-geog.dist[clusters$kelpbio==1,clusters$kelpbio==1]
geog.dist.c2<-geog.dist[clusters$kelpbio==2,clusters$kelpbio==2]

cormat.lb.c1<-cormat.lb.filt[clusters$kelpbio==1,clusters$kelpbio==1]
cormat.ub.c1<-cormat.ub.filt[clusters$kelpbio==1,clusters$kelpbio==1]
cormat.lb.c2<-cormat.lb.filt[clusters$kelpbio==2,clusters$kelpbio==2]
cormat.ub.c2<-cormat.ub.filt[clusters$kelpbio==2,clusters$kelpbio==2]

sncf.lb.c1<-splineFit(geog.dist[clusters$kelpbio==1,clusters$kelpbio==1], cormat.lb.c1, quantiles=c(0.025,0.5,0.975))
sncf.ub.c1<-splineFit(geog.dist[clusters$kelpbio==1,clusters$kelpbio==1], cormat.ub.c1, quantiles=c(0.025,0.5,0.975))
sncf.lb.c2<-splineFit(geog.dist[clusters$kelpbio==2,clusters$kelpbio==2], cormat.lb.c2, quantiles=c(0.025,0.5,0.975))
sncf.ub.c2<-splineFit(geog.dist[clusters$kelpbio==2,clusters$kelpbio==2], cormat.ub.c2, quantiles=c(0.025,0.5,0.975))


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

pdf("./TailAssociation/dist_decay_central.pdf")
plot(sncf.lb.c1$emp.spline$x, sncf.lb.c1$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
     xlab="Distance (km)", ylab="Partial Spearman correlation", main="Central California (c1)")
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[1,], col="red", lty=3)
# lines(sncf.lb.c1$emp.spline$x, sncf.lb.c1$spline.quantiles[3,], col="red", lty=3)
lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[1,], col="blue", lty=3)
# lines(sncf.ub.c1$emp.spline$x, sncf.ub.c1$spline.quantiles[3,], col="blue", lty=3)
legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
dev.off()

pdf("./TailAssociation/dist_decay_southern.pdf")
plot(sncf.lb.c2$emp.spline$x, sncf.lb.c2$emp.spline$y, col="red", type="l", ylim=c(0,0.4), lwd=2,
     xlab="Distance (km)", ylab="Partial Spearman correlation", main="Southern California (c2)")
# lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[1,], col="red", lty=3)
# lines(sncf.lb.c2$emp.spline$x, sncf.lb.c2$spline.quantiles[3,], col="red", lty=3)
lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$emp.spline$y, col="blue", lwd=2)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[1,], col="blue", lty=3)
# lines(sncf.ub.c2$emp.spline$x, sncf.ub.c2$spline.quantiles[3,], col="blue", lty=3)
legend("topright",lty=1,col=c("red","blue"),legend=c("lower","upper"))
dev.off()


## visualize differences between upper and lower tails---------------------------------------------

cormat.diff<-cormat.lb.filt-cormat.ub.filt

pdf("./TailAssociation/hist_kelpsynch_taildiff.pdf")
hist(cormat.diff[lower.tri(cormat.diff)])
abline(v=median(cormat.diff[lower.tri(cormat.diff)], na.rm=T), lty=2, col="red")
dev.off()

mean(cormat.diff[lower.tri(cormat.diff)])


cormat.diff.nodiag<-cormat.diff
diag(cormat.diff.nodiag)<-NA
pdf("./TailAssociation/synchmat_kelp_taildiff.pdf")
image.plot(cormat.diff.nodiag)
dev.off()



## test whether geographies of kelp synchrony differ between upper and lower tails ----------------
mantel(as.dist(cormat.lb)~as.dist(cormat.ub))


# image(cormat.lb)
# image(cormat.ub)


matrixDiff<-function(mat1, mat2, tslength, nreps=100){
  library(mvtnorm)
  
  emp.diff <- 1-cor(mat1[lower.tri(mat1)],mat2[lower.tri(mat2)], use="pairwise.complete.obs")

  covmat<-mat1  
  diag(covmat)<-1
  
  surrog.diff <- rep(NA, nreps)
  for(rep in 1:nreps){
    surrdat<-rmvnorm(n=tslength, sigma = covmat)
    surrcor<-cor(surrdat)
    surrog.diff[rep]<-1-cor(mat1[lower.tri(mat1)],surrcor[lower.tri(surrcor)], use="pairwise.complete.obs")

  }
  
  surrog.distn<-ecdf(surrog.diff)
  p<-1-surrog.distn(emp.diff)
  
  return(list(emp.diff=emp.diff, p_value=p))
}

matrixDiff(cormat.ub, cormat.lb, 33, 1000)



## Measure average synchrony within a distance threshold

avg_by_dist<-function(cormat, dmat, thresh){
  
  diag(cormat)<-NA
  
  nbhd<-ifelse(dmat <= thresh, 1, 0)
  nbhd.cor<-cormat * nbhd
  return(apply(nbhd.cor, 1, mean, na.rm=T))
  
}

geog.dist<-ncf::gcdist(clusters$Lon, clusters$Lat)

coravg.lb.50km<-avg_by_dist(cormat.lb, geog.dist, 50)
coravg.ub.50km<-avg_by_dist(cormat.ub, geog.dist, 50)

coravg.diff.50km<-abs(coravg.lb.50km)-abs(coravg.ub.50km)

nclasses=7
pal<-brewer.pal(nclasses,"RdYlBu")

class_eq_int<-function(x, n){
  
  breaks<-seq(from=min(x, na.rm=T), to=max(x, na.rm=T)+1/length(x), length.out=n+1)
  class<-rep(NA, length(x))
  
  for(ii in 1:n){
    class[x>=breaks[ii] & x<breaks[ii+1]]<-ii
  }
  return(list(class=class, breaks=breaks))
}


pdf("./TailAssociation/sync_taildiff_50km.pdf")
plot(clusters$Lon, clusters$Lat, pch=19, col=pal[class_eq_int(coravg.diff.50km,nclasses)$class], cex=0.5,
     xlab="Longitude", ylab="Latitude", main="mean <=50 km synchrony lower tail - upper tail")
legend("topright", pch=19, legend=c("-0.02 to -0.01","-0.01 to 0.00","0.00 to 0.02","0.02 to 0.03","0.03 to 0.04","0.04 to 0.05","0.05 to 0.07"),
       col=pal)
dev.off()


#Tail dependence in environmental relationships

kelpXwaves.all<-rep(NA, nrow(kelp.ann))
kelpXwaves.lb<-rep(NA, nrow(kelp.ann))
kelpXwaves.ub<-rep(NA, nrow(kelp.ann))
kelpXno3.all<-rep(NA, nrow(kelp.ann))
kelpXno3.lb<-rep(NA, nrow(kelp.ann))
kelpXno3.ub<-rep(NA, nrow(kelp.ann))
kelpXnpgo.lb<-rep(NA, nrow(kelp.ann))
kelpXnpgo.ub<-rep(NA, nrow(kelp.ann))
kelpXnpgo.all<-rep(NA, nrow(kelp.ann))

for(ii in 1:nrow(kelp.ann)){
  
  kelpXwaves.lb[ii]<-partialSpearman(kelp.ann[ii,],-1*waves.ann[ii,],c(0,0.5))
  kelpXwaves.ub[ii]<-partialSpearman(kelp.ann[ii,],-1*waves.ann[ii,],c(0.5,1))
  kelpXwaves.all[ii]<-cor(kelp.ann[ii,],-1*waves.ann[ii,], method="spearman")
  kelpXno3.lb[ii]<-partialSpearman(kelp.ann[ii,],no3.ann[ii,],c(0,0.5))
  kelpXno3.ub[ii]<-partialSpearman(kelp.ann[ii,],no3.ann[ii,],c(0.5,1))
  kelpXno3.all[ii]<-cor(kelp.ann[ii,],no3.ann[ii,], method="spearman")
  kelpXnpgo.lb[ii]<-partialSpearman(kelp.ann[ii,],npgo.ann,c(0,0.5))
  kelpXnpgo.ub[ii]<-partialSpearman(kelp.ann[ii,],npgo.ann,c(0.5,1))
  kelpXnpgo.all[ii]<-cor(kelp.ann[ii,],no3.ann[ii,], method="spearman")
  
}


kelpXwaves.filt<-ifelse(kelpXwaves.all>=0,1,NA)
kelpXno3.filt<-ifelse(kelpXno3.all>=0,1,NA)
kelpXnpgo.filt<-ifelse(kelpXnpgo.all>=0,1,NA)

kelpXwaves.diff<-kelpXwaves.lb-kelpXwaves.ub*kelpXwaves.filt
kelpXno3.diff<-kelpXno3.lb-kelpXno3.ub*kelpXno3.filt
kelpXnpgo.diff<-kelpXnpgo.lb-kelpXnpgo.ub*kelpXnpgo.filt

pdf("./TailAssociation/hist_taildiff_kelpXwaves.pdf")
hist(kelpXwaves.diff, xlab="lower tail - upper tail")
abline(v=median(kelpXwaves.diff, na.rm=T), col="red", lty=2)
dev.off()

pdf("./TailAssociation/hist_taildiff_kelpXno3.pdf")
hist(kelpXno3.diff, xlab="lower tail - upper tail")
abline(v=median(kelpXno3.diff,na.rm=T), col="red", lty=2)
dev.off()

pdf("./TailAssociation/hist_taildiff_kelpXnpgo.pdf")
hist(kelpXnpgo.diff, xlab="lower tail - upper tail")
abline(v=median(kelpXnpgo.diff,na.rm=T), col="red", lty=2)
dev.off()

cor(kelpXwaves.diff, kelpXno3.diff, use="pairwise.complete.obs") #these are slightly positively correlated

cor(kelpXwaves.diff, kelpXnpgo.diff, use="pairwise.complete.obs") #these are moderately positively correlated

cor(kelpXno3.diff, kelpXnpgo.diff, use="pairwise.complete.obs") #these are moderately positively correlated

pdf("./TailAssociation/map_taildiff_kelpXwaves.pdf")
plot(clusters$Lon, clusters$Lat, pch=19, col=pal[class_eq_int(kelpXwaves.diff,nclasses)$class], cex=0.5,
     xlab="Longitude", ylab="Latitude", main="Waves lower tail - upper tail")
legend("topright", pch=19, legend=c("-0.27 to -0.19","-0.19 to -0.10","-0.10 to -0.02","-0.02 to 0.07","0.07 to 0.15",
                                    "0.15 to 0.23","0.23 to 0.32"),
       col=pal)
dev.off()

pdf("./TailAssociation/map_taildiff_kelpXno3.pdf")
plot(clusters$Lon, clusters$Lat, pch=19, col=pal[class_eq_int(kelpXno3.diff,nclasses)$class], cex=0.5,
     xlab="Longitude", ylab="Latitude", main="NO3 lower tail - upper tail")
legend("topright", pch=19, legend=c("-0.20 to -0.14","-0.14 to -0.08","-0.08 to -0.02","-0.02 to 0.04","0.04 to 0.10",
                                    "0.10 to 0.16","0.16 to 0.22"),
       col=pal)
dev.off()

pdf("./TailAssociation/map_taildiff_kelpXnpgo.pdf")
plot(clusters$Lon, clusters$Lat, pch=19, col=pal[class_eq_int(kelpXnpgo.diff,nclasses)$class], cex=0.5,
     xlab="Longitude", ylab="Latitude", main="NPGO lower tail - upper tail")
legend("topright", pch=19, legend=c("-0.30 to -0.23","-0.23 to -0.17","-0.17 to -0.10","-0.10 to -0.03","-0.03 to 0.04",
                                    "0.04 to 0.10","0.10 to 0.17"),
       col=pal)
dev.off()



#make the environmental covariates averages within 50km!

avg_by_dist2<-function(x, dmat, thresh){
  
  nbhd<-ifelse(dmat <= thresh, 1, 0)
  nbhd.cor<-x * nbhd
  return(apply(nbhd.cor, 1, mean, na.rm=T))
  
}




moddat<-data.frame(kelp = coravg.diff.50km,
                   no3 = avg_by_dist2(kelpXno3.diff, geog.dist, 50),
                   waves = avg_by_dist2(kelpXwaves.diff, geog.dist, 50),
                   npgo = avg_by_dist2(kelpXnpgo.diff, geog.dist, 50),
                   lat = clusters$Lat,
                   lon = clusters$Lon)

moddat<-moddat[complete.cases(moddat),]

mod<-gls(kelp ~ no3 + waves + npgo, correlation = corExp(form = ~ lon + lat), data=moddat)
summary(mod)
