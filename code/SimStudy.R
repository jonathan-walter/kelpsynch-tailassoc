#This script carries out the simulation study described in the sup mat text that Dan wrote called "General theory".
#Reuman
#2021 07 19

#***Preliminaries

#Libraries needed to run. Namespaces not loaded, since :: is always used, but these need to be installed.
#testthat - needed for the embedded testing only
#copula
#parallel

#Cores on your machine you can spare. This code uses mclapply, works on linus, does not really work on
#widows, seems likely to work on mac
cores_to_use<-10 #number of cores you can spare

#***Functions

#Function for simulating the model for given parameters
#
#Args
#A          Autoregressive parameter, as described in the sup mat section that Dan wrote called "General theory".
#sigma_l    Standard deviation parameter for local noise. See "General theory" section.
#delta_bar  Mean parameter for synchronized noise. See "General theory" section.
#sigma_s    Standard deviation parameter for synchronized noise. See "General theory" section.
#burnin     Length of burnin in time steps
#simlen     Number of time steps to keep after the burnin
#
#Output
#simout     Output of the simulations, 2 by simlen matrix
#
runsim<-function(A,sigma_l,delta_bar,sigma_s,burnin,simlen)
{
  #error checking parameters
  if (A<(-1) || A>1)
  {
    stop("Error in runsim: parameter A out of range")
  }
  if (sigma_l<=0)
  {
    stop("Error in runsim: parameter sigma_l out of range")
  }
  if (sigma_s<=0)
  {
    stop("Error in runsim: parameter sigma_s out of range")
  }
  
  #sim length
  totsimlen<-burnin+simlen
  
  #get the necessary synchronous noise
  syncnoise<-rnorm(totsimlen-1,delta_bar,sigma_s)
  syncnoise<-1/(1+exp(-syncnoise))
  syncnoise<-syncnoise-mean(syncnoise)
  
  #get the necessary local noise
  locnoise<-matrix(rnorm(2*(totsimlen-1),0,sigma_l),2,totsimlen-1)
  
  #now do the sim
  res<-matrix(NA,2,totsimlen)
  res[,1]<-0
  for (counter in 2:totsimlen)
  {
    res[,counter]<-A*res[,counter-1]+syncnoise[counter-1]+locnoise[,counter-1]
  }
  
  #cut the burnin and return
  res<-res[,(burnin+1):totsimlen]
  return(res)
}

#Some very modest testing of the sim function
#res<-runsim(A=.5,sigma_l=1,delta_bar=1,sigma_s=1,burnin=100,simlen=1000)
#testthat::expect_equal(class(res),"matrix")
#testthat::expect_equal(dim(res),c(2,1000))
#testthat::expect_equal(sum(is.na(res)),0)

#Function for partial Spearman correlation
#
#Args
#x, y         Numeric vectors of the same length between which the partial Spearman correlation is computed.
#               Assumed that marginals have already been normalized, if necessary (e.g., used copula::pobs).
#lb, ub       Lower and upper bounds between 0 and 1
#
#Output
#The partial Spearman correlation
#
#Notes: For speed, no error checking done. For instance, if you pass the original variables instead of their
#normalized ranks, it'll give you the wrong answer with no warning.
#
partial_spearman<-function(x,y,lb,ub)
{
  #get mean and variance
  x_mean<-mean(x)
  y_mean<-mean(y)
  x_var<-var(x)
  y_var<-var(y)
  
  #compute the indices of the points between the bounds
  inds<-which(x+y>2*lb & x+y<2*ub)
  
  #get the portion of the Spearman
  res<-sum((x[inds]-x_mean)*(y[inds]-y_mean))/((length(x)-1)*sqrt(x_var*y_var))
  
  return(res)  
}

#some modest testing of the partial_spearman function
#ccop<-copula::claytonCopula(2,2)
#d<-copula::rCopula(1000,ccop)
#d<-copula::pobs(d)
#rl<-partial_spearman(d[,1],d[,2],0,.5)
#rh<-partial_spearman(d[,1],d[,2],.5,1) 
#expect_true(rl>rh) #the lower-tail value should be substantially bigger
#rl
#rh
#ra<-cor(d[,1],d[,2],method="spearman")
#rl+rh
#ra
#expect_equal(ra,rl+rh) #note only works because we used pobs on the data, d, above

#***Now do the simulation study

#the values we will use for each parameter, all combinations of values to be used
all_A<-c(-.8,-.5,-.3,0,.3,.5,.8)
all_sigma_l<-c(.1,.2,.3)
all_delta_bar<-seq(from=-8,to=8,by=.1)
all_sigma_s<-seq(from=.5,to=3,by=.1)
burnin<-100
simlen<-1000

#make an mclapply arglist, to be passed as the X argument of mclapply
mclapply_arglist<-list()
len<-0
for (A_counter in 1:length(all_A))
{
  for (sigma_l_counter in 1:length(all_sigma_l))
  {
    h<-list(A=all_A[A_counter],
            sigma_l=all_sigma_l[sigma_l_counter],
            all_delta_bar=all_delta_bar,
            all_sigma_s=all_sigma_s,
            burnin=burnin,
            simlen=simlen)
    mclapply_arglist[[len+1]]<-h
    len<-len+1
  }
}

#make wrapper function to be passed as the FUN argument of mclapply
mclapply_wrapper<-function(x)
{
  #unpack arguments
  A<-x$A
  sigma_l<-x$sigma_l
  all_delta_bar<-x$all_delta_bar
  all_sigma_s<-x$all_sigma_s
  burnin<-x$burnin
  simlen<-x$simlen
  
  wrapper_res<-matrix(NA,length(all_delta_bar),length(all_sigma_s))
  for (delta_bar_counter in 1:length(all_delta_bar))
  {
    for (sigma_s_counter in 1:length(all_sigma_s))
    {
      simres<-runsim(A,sigma_l,
                     all_delta_bar[delta_bar_counter],
                     all_sigma_s[sigma_s_counter],
                     burnin=burnin,simlen=simlen)
      simres<-copula::pobs(t(simres))
      wrapper_res[delta_bar_counter,sigma_s_counter]<-partial_spearman(simres[,1],simres[,2],.5,1)-partial_spearman(simres[,1],simres[,2],0,.5)
    }
  }
    
  return(wrapper_res)
}

#run it once to test, and to time it
#system.time(h<-mclapply_wrapper(mclapply_arglist[[1]]))
#class(h)
#dim(h)
#length(all_delta_bar)
#length(all_sigma_s)
#sum(is.na(h))
#range(h)

#now do the while sim study
ssres<-parallel::mclapply(mclapply_arglist,mclapply_wrapper,mc.cores=cores_to_use)

#***Now plot results

#get the limits of the values
zlimits<-c(0,0)
for (counter in 1:length(ssres))
{
  zlimits<-range(zlimits,ssres[[counter]])
}
h<-max(abs(zlimits))
zlimits<-c(-h,h) #to make the symmetric about 0

#technical stuff for the split colorbar
breaks<-seq(from=zlimits[1],to=zlimits[2],length.out=102)
colbartop<-rgb(red=1,green=seq(from=1,to=0,length.out=51),blue=seq(from=1,to=0,length.out=51))
#plot(1:51,1:51,col=colbartop,pch=20,cex=2)
colbarbot<-rgb(blue=1,green=seq(from=0,to=1,length.out=51),red=seq(from=0,to=1,length.out=51))
#plot(1:51,1:51,col=colbarbot,pch=20,cex=2)
colbar<-c(colbarbot[1:50],colbartop)
#plot(1:101,1:101,col=colbar,pch=20,cex=2)

#plotting dimensions, in inches
panht<-1
panwd<-panht
xaxht<-.8
yaxwd<-.8
gap<-.1
leght<-.5
totht<-xaxht+length(all_A)*(panht+gap)+leght+gap
totwd<-yaxwd+length(all_sigma_l)*(panwd+gap)
pdf(file="SimStudyFig.pdf",width=totwd,height=totht)

len<-1
for (A_counter in 1:length(all_A))
{
  for (sigma_l_counter in 1:length(all_sigma_l))
  {
    #set up the panel for these values of A and sigma_l
    if (A_counter==1 && sigma_l_counter==1)
    {
      par(fig=c((yaxwd+(sigma_l_counter-1)*(panwd+gap))/totwd,
                (yaxwd+(sigma_l_counter-1)*(panwd+gap)+panwd)/totwd,
                (xaxht+(A_counter-1)*(panht+gap))/totht,
                (xaxht+(A_counter-1)*(panht+gap)+panht)/totht),
          mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25)
    } else
    {
      par(fig=c((yaxwd+(sigma_l_counter-1)*(panwd+gap))/totwd,
                (yaxwd+(sigma_l_counter-1)*(panwd+gap)+panwd)/totwd,
                (xaxht+(A_counter-1)*(panht+gap))/totht,
                (xaxht+(A_counter-1)*(panht+gap)+panht)/totht),
          mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
    }
    
    #make the plot
    image(all_delta_bar,all_sigma_s,ssres[[len]],breaks=breaks,col=colbar,xaxt="n",yaxt="n")
    len<-len+1
    lines(c(0,0),range(all_sigma_s))
    
    #create the axes
    if (A_counter==1)
    {
      axis(side=1)
      mtext(expression(bar(italic(delta))),side=1,line=1.2)
      #mtext(expression(paste(italic(sigma)[italic(l)],"=")),side=1,line=2.9)
      mtext(bquote(italic(sigma)[italic(l)]==.(all_sigma_l[sigma_l_counter])),side=1,line=2.9)
    }
    if (sigma_l_counter==1)
    {
      axis(side=2)
      mtext(expression(italic(sigma)[italic(s)]),side=2,line=1.2)
      mtext(bquote(italic(A)==.(all_A[A_counter])),side=2,line=2.9)
    }
  }
}

#now make the legend
par(fig=c((gap)/totwd,
          (totwd-gap)/totwd,
          (xaxht+length(all_A)*(panht+gap))/totht,
          (xaxht+length(all_A)*(panht+gap)+leght)/totht),
    mai=c(0,0,0,0),mgp=c(3,.15,0),tcl=-.25,new=TRUE)
plot(c(0,1),c(0,1),lty="n")
numbox<-length(colbar)
rect(xleft=seq(from=0,to=1-1/numbox,length.out=numbox),
     xright=seq(from=1/numbox,to=1,length.out=numbox),
     ybottom=2/3,
     ytop=1,
     col=colbar,border=NA)
ticlocs<-c(-.25,-.15,0,.15,.25)
ep<-.1
for (counter in 1:length(ticlocs))
{
  lines(rep((ticlocs[counter]-zlimits[1])/diff(zlimits),2),c(2/3,2/3-ep),type="l")
  text((ticlocs[counter]-zlimits[1])/diff(zlimits),2/3-3*ep,as.character(ticlocs[counter]))
}

dev.off()
