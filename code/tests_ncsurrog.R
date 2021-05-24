rm(list=ls())

library(testthat)
source("ncsurrog.R")

theseed<-101

#test the error catching
set.seed(theseed)
m<-matrix(rnorm(1000),100,10)
m[1,1]<-NA
numsurrog<-10
plotcheckloc<-NA
expect_error(ncsurrog(m=m,numsurrog=numsurrog,plotcheckloc=plotcheckloc),"Error in ncsurrog: NA, Inf, NaN, etc. not allowed in m")

#test the correct format of output
set.seed(theseed)
m<-matrix(rnorm(1000),100,10)
numsurrog<-15
plotcheckloc<-NA
surm<-ncsurrog(m=m,numsurrog=numsurrog,plotcheckloc=plotcheckloc)
expect_equal(class(surm),"array")
expect_equal(dim(surm),c(100,10,numsurrog))

#check that, for each surrogate and each location (OK, I only did a few), the surrogate
#values for that location are a rearrangement of the actual data values
expect_equal(sort(m[,1]),sort(surm[,1,1]))
expect_equal(sort(m[,1]),sort(surm[,1,2]))
expect_equal(sort(m[,1]),sort(surm[,1,3]))
expect_equal(sort(m[,2]),sort(surm[,2,1]))
expect_equal(sort(m[,2]),sort(surm[,2,2]))
expect_equal(sort(m[,2]),sort(surm[,2,3]))
expect_equal(sort(m[,3]),sort(surm[,3,1]))
expect_equal(sort(m[,3]),sort(surm[,3,2]))
expect_equal(sort(m[,3]),sort(surm[,3,3]))

#***Some of the following checks are visual, not automated, checking for reasonableness of results

#make some data with asymmetric tail association
set.seed(theseed)
ccop<-claytonCopula(param=3,dim=10)
dat<-rCopula(1000,ccop)
plot(dat[,1],dat[,2],type="p")
plot(dat[,3],dat[,4],type="p")

#get surrogates and check for a few of them at a few locations that surrogates are premutations of the data
surm<-ncsurrog(m=dat,numsurrog=1000,plotcheckloc="./test")
expect_equal(sort(surm[,1,1]),sort(dat[,1]))
expect_equal(sort(surm[,1,2]),sort(dat[,1]))
expect_equal(sort(surm[,1,3]),sort(dat[,1]))
expect_equal(sort(surm[,2,1]),sort(dat[,2]))
expect_equal(sort(surm[,2,2]),sort(dat[,2]))
expect_equal(sort(surm[,2,3]),sort(dat[,2]))
expect_equal(sort(surm[,3,1]),sort(dat[,3]))
expect_equal(sort(surm[,3,2]),sort(dat[,3]))
expect_equal(sort(surm[,3,3]),sort(dat[,3]))

#check the surrogates look like they have symmetric tail association, in contrast to the original data
plot(dat[,1],dat[,2],type="p")
plot(surm[,1,1],surm[,2,1],type="p")
plot(surm[,1,2],surm[,2,2],type="p")
plot(surm[,1,3],surm[,2,3],type="p")
plot(surm[,1,4],surm[,2,4],type="p")
plot(surm[,1,5],surm[,2,5],type="p")
plot(surm[,1,6],surm[,2,6],type="p")

plot(surm[,3,1],surm[,6,1],type="p")
plot(surm[,3,2],surm[,6,2],type="p")
plot(surm[,3,3],surm[,6,3],type="p")
plot(surm[,3,4],surm[,6,4],type="p")
plot(surm[,3,5],surm[,6,5],type="p")
plot(surm[,3,6],surm[,6,6],type="p")

#you can also open the plotcheck result saved in ./test.pdf and make sure the dashed line passes 
#through the quantiles.