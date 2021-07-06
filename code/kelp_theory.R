## theoretical example for kelp manuscript

rm(list=ls())

burn = 100
tmax = burn+50
sd.reg = 1
mu.reg1 = -2 # mean of regional driver for location 1
mu.reg2 = 2 # mean of regional driver for location 2
sd.loc = 0.1 # sd of local driver (mean = 0)
n.each = 2
r = 0.9
K = 100

set.seed(19)

#sigmoid function acts as filter on environmental noises
sigmoid <- function(x) {
  return(1/(1+exp(-x)))
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

Nt.1[,1] <- rnorm(n.each, mean = K, sd = K/10)
Nt.2[,1] <- rnorm(n.each, mean = K, sd = K/10)


for(tt in 2:tmax){
  
  Nt.1[,tt] <- ricker(Nt.1[,tt-1], r, K) + rescale(sigmoid(Et.1), K/10)[tt] + rnorm(n.each, sd=K/20)
  Nt.2[,tt] <- ricker(Nt.2[,tt-1], r, K) + rescale(sigmoid(Et.2), K/10)[tt] + rnorm(n.each, sd=K/20)
  
}

plot(Nt.1[1,(burn+1):tmax], Nt.1[2,(burn+1):tmax])
plot(Nt.2[1,(burn+1):tmax], Nt.2[2,(burn+1):tmax])


plot(rep(Et.1[(burn+1):tmax], each=2), Nt.1[,(burn+1):tmax])
plot(rep(Et.2[(burn+1):tmax], each=2), Nt.2[,(burn+1):tmax])


plot(Nt.1[1, (burn+1):tmax], type="l")
lines(Nt.1[2,(burn+1):tmax], lty=2)


plot(Nt.2[1, (burn+1):tmax], type="l")
lines(Nt.2[2,(burn+1):tmax], lty=2)
