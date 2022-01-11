library(copula)
library(mvtnorm)

set.seed(11)

rho <- 0.9
np <- 500

# for clayton : left tail dep.
copC<-claytonCopula(4)
parC<-iRho(copC,rho=rho)

# first generate Clayton copula 
cc<-claytonCopula(par=parC,dim=2)
mtc_l<-rCopula(np,cc)  

# flip Clayton copula 
mtc_r<-1-mtc_l

noATA <- rmvnorm(np, sigma = matrix(c(1,rho,rho,1),2,2))


png("~/GitHub/kelpsynch-tailassoc/manuscript/fig1_pedagogical.png", units="in",
    res=300, width=6.5, height=2.5)

par(mfrow=c(1,3), mar=c(3.1,3.1,1.6,1.1), mgp=c(1.5,1,0))

plot(qnorm(mtc_l[,1]), qnorm(mtc_l[,2]), xaxt="n", yaxt="n",
     xlab="Population 1", ylab="Population 2", pch=16, cex=0.7)
axis(1, labels=NA); axis(2, labels=NA)
mtext("Lower tail dependence", line=0.25, cex=0.8)
text(par("usr")[1] + 0.05*diff(par("usr")[1:2]), 
     par("usr")[4] - 0.05*diff(par("usr")[3:4]), "a)")

plot(noATA[,1], noATA[,2], xaxt="n", yaxt="n",
     xlab="Population 1", ylab="Population 2", pch=16, cex=0.7)
axis(1, labels=NA); axis(2, labels=NA)
mtext("No tail dependence", line=0.25, cex=0.8)
text(par("usr")[1] + 0.05*diff(par("usr")[1:2]), 
     par("usr")[4] - 0.05*diff(par("usr")[3:4]), "b)")

plot(qnorm(mtc_r[,1]), qnorm(mtc_r[,2]), xaxt="n", yaxt="n",
     xlab="Population 1", ylab="Population 2", pch=16, cex=0.7)
axis(1, labels=NA); axis(2, labels=NA)
mtext("Upper tail dependence", line=0.25, cex=0.8)
text(par("usr")[1] + 0.05*diff(par("usr")[1:2]), 
     par("usr")[4] - 0.05*diff(par("usr")[3:4]), "c)")

dev.off()