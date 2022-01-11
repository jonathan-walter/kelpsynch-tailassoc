## theoretical example for kelp manuscript

rm(list=ls())

burn = 50
tmax = burn+100
sd.reg = 1.5
mu.reg1 = -2 # mean of regional driver for location 1
mu.reg2 = 2 # mean of regional driver for location 2
sd.loc = 0.15 # sd of local driver (mean = 0)
n.each = 2
b = 0.5
c=1.5

set.seed(17)

#sigmoid function acts as filter on environmental noises
sigmoid <- function(x,c=1) {
  return(1/(1+exp(-c*x)))
}

x = seq(-6,6,0.2)
plot(x, sigmoid(x))

#simple ricker population growth
ricker <- function(N, r, K){ 
  return(N*exp(r*(1-(N/K))))
}

#set a time series to have mean zero and specified sd
rescale <- function(x, newsd){
  y <- x-mean(x)
  y <- y/(sd(x)/newsd)
  return(y)
}




Nt.1 <- matrix(NA, nrow=n.each, ncol=tmax)
Nt.2 <- Nt.1

Et.1 <- rnorm(tmax, mean=mu.reg1, sd=sd.reg)
Et.2 <- rnorm(tmax, mean=mu.reg2, sd=sd.reg)

Nt.1[,1] <- rnorm(n.each, mean = 0, sd = 1)
Nt.2[,1] <- rnorm(n.each, mean = 0, sd = 1)


for(tt in 2:tmax){
  
  Nt.1[,tt] <- Nt.1[,tt-1]*b + (sigmoid(Et.1, c)[tt] - sigmoid(mu.reg1, c)) + rnorm(n.each, sd=sd.loc)
  Nt.2[,tt] <- Nt.2[,tt-1]*b + (sigmoid(Et.2, c)[tt] - sigmoid(mu.reg2, c)) + rnorm(n.each, sd=sd.loc)
  
}

plot(Nt.1[1,(burn+1):tmax], Nt.1[2,(burn+1):tmax])
plot(Nt.2[1,(burn+1):tmax], Nt.2[2,(burn+1):tmax])


plot(rep(Et.1[(burn+1):tmax], each=2), Nt.1[,(burn+1):tmax])
plot(rep(Et.2[(burn+1):tmax], each=2), Nt.2[,(burn+1):tmax])


plot(Nt.1[1, (burn+1):tmax], type="l")
lines(Nt.1[2,(burn+1):tmax], lty=2)


plot(Nt.2[1, (burn+1):tmax], type="l")
lines(Nt.2[2,(burn+1):tmax], lty=2)





laymat <- matrix(NA, nrow=3, ncol=8)
laymat[1:2,1:4] <- 1
laymat[1,5:8] <- c(2,2,3,3)
laymat[2,5:8] <- c(4,4,5,5)
laymat[3,] <- rep(c(6,7),each=4)


png("~/GitHub/kelpsynch-tailassoc/manuscript/fig2_theory_ar1.png", units="in", res=300, width=6.5, height=4.5)

layout(laymat)
par(mar=c(2.3,2.3,1.1,1.1), mgp=c(1,0.5,0))

#first panel
plot(x, sigmoid(x,c), type="l", lwd=2, xlab = "Environmental driver", ylab="Population", ylim=c(-0.1,1.1),
     xaxt="n", yaxt="n")
axis(1, at=seq(min(x),max(x),length.out=5), labels=FALSE)
axis(2, at=seq(0,1,length.out=5),labels=FALSE)
arrows(x0=mu.reg1-1.96*sd.reg, y0=-0.05, x1=mu.reg1+1.96*sd.reg, y1=-0.05, col="blue", length=0.04, angle=90, code=3, lwd=1.3)
text(mu.reg1,-0.11,"Effects stronger in upper tail",col="blue")
text(mu.reg2,1.1,"Effects stronger in lower tail",col="red")
arrows(x0=mu.reg2-1.96*sd.reg, y0=1.05, x1=mu.reg2+1.96*sd.reg, y1=1.05, col="red", length=0.04, angle=90, code=3, lwd=1.3)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.05*abs(diff(par("usr")[3:4])),"a)")

#second panel
plot(rep(Et.1[(burn+1):tmax], each=2), Nt.1[,(burn+1):tmax], pch=rep(c(16,1),times=length((burn+1):tmax)),
     col="blue", xlab="Environmental driver", ylab="Population", xaxt="n", yaxt="n")
axis(1, labels=FALSE)
axis(2, labels=FALSE)
text(par("usr")[1]+0.1*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"b)")
mtext("Upper tail dependence", cex=0.6, line=0.1)
legend("bottomright", pch=c(16,1), legend=c("Pop 1", "Pop 2"), bty="n", cex=0.8)

#third panel
plot(rep(Et.2[(burn+1):tmax], each=2), Nt.2[,(burn+1):tmax],  pch=rep(c(16,1),times=length((burn+1):tmax)),
     col="red", xlab="Environmental driver", ylab="Population", xaxt="n", yaxt="n")
axis(1, labels=FALSE)
axis(2, labels=FALSE)
text(par("usr")[1]+0.1*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"c)")
mtext("Lower tail dependence", cex=0.6, line=0.1)

#fourth panel
plot(Nt.1[1,(burn+1):tmax], Nt.1[2,(burn+1):tmax], pch=16, col="blue", xlab="Population 1",
     ylab="Population 2", xaxt="n", yaxt="n")
axis(1, labels=FALSE)
axis(2, labels=FALSE)
text(par("usr")[1]+0.1*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"d)")

#fifth panel
plot(Nt.2[1,(burn+1):tmax], Nt.2[2,(burn+1):tmax], pch=16, col="red", xlab="Population 1",
     ylab="Population 2", xaxt="n", yaxt="n")
axis(1, labels=FALSE)
axis(2, labels=FALSE)
text(par("usr")[1]+0.1*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"e)")

#sixth panel
plot(Nt.1[1, (burn+1):tmax], type="l", col="blue", xlab="Time", ylab="Population", xaxt="n", yaxt="n",
     #ylim=c(min(Nt.1[,(burn+1):tmax]), max(Nt.1[,(burn+1):tmax])))
     ylim=c(-0.7,1.5))
axis(1, labels=FALSE)
axis(2, labels=FALSE)
lines(Nt.1[2,(burn+1):tmax], lty=2, col="blue")
segments(x0=37.5, y0=1.35, x1=c(3,11), y1=c(1.22, 1.33), col="darkgrey")
segments(x0=37.5, y0=-0.5, x1=25, y1=-0.46, col="darkgrey")
segments(x0=76, y0=-0.5, x1=88, y1=-0.37, col="darkgrey")
text(37.5, 1.35, "Synchronous booms", pos=4, col="darkgrey", cex=0.8)
text(37.5, -0.5, "Asynchronous crashes", pos=4, col="darkgrey", cex=0.8)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"f)")

#seventh
plot(Nt.2[1, (burn+1):tmax], type="l", col="red", xlab="Time", ylab="Population", xaxt="n", yaxt="n",
     #ylim=c(min(Nt.2[,(burn+1):tmax]), max(Nt.2[,(burn+1):tmax])))
     ylim=c(-1.2,0.7))
axis(1, labels=FALSE)
axis(2, labels=FALSE)
lines(Nt.2[2,(burn+1):tmax], lty=2, col="red")
segments(x0=16, y0=0.65, x1=6, y1=0.59, col="darkgrey")
segments(x0=54, y0=0.65, x1=65, y1=0.48, col="darkgrey")
segments(x0=50, y0=-1.1, x1=36, y1=-1.05, col="darkgrey")
segments(x0=63, y0=-0.95, x1=73, y1=-0.8, col="darkgrey")
text(16,0.65, "Asynchronous booms", pos=4, col="darkgrey", cex=0.8)
text(50,-1.1, "Synchronous crashes", pos=4, col="darkgrey", cex=0.8)
text(par("usr")[1]+0.05*abs(diff(par("usr")[1:2])), par("usr")[4]-0.1*abs(diff(par("usr")[3:4])),"g)")


dev.off()


