#Packages which must be installed to run this: 
#copula
#mvtnorm
#matrixcalc
#Matrix
#grDevices

#Another function which is a needed:
source("./code/alignranks.R")

#This function takes a bunch of time series measured at the same
#times in different locations and creates surrogate datasets which
#are statistically similar except the copula has been randomized
#to a copula in the normal family. "Statistcally similar" means the 
#location marginals are (exactly) the same and the pairwise spearman 
#correlations of time series between locations are similar (differing
#only by a kind of sampling variation).
#
#Notes: In some cases the algorithm does not strictly work because an 
#intermediate matrix which needs to be positive semi-definite is not. 
#In that case a warning is given and a "nearby" positive semi-definite
#matrix is used instead. This is risk and requires some due dilligence,
#so speak to Dan if you get this warning. There is a line of code you
#can uncomment if you just want to throw an error here, instead.
#
#Args
#m              A N by n matrix, where N is the length of the time 
#                 series and n is the number of time series
#numsurrog      The desired number of surrogate datasets
#plotcheckloc   Location to store a pdf giving a visual check of how 
#                 close surrogate Spearman values are to real data 
#                 Spearman values. Default NA means skip it. The 
#                 filename can be a path, and should not have the 
#                 extension ".pdf". The plot shows real-data Spearman
#                 values for each pair of locations on the x axis, 
#                 and the corresponging distribution of surrogate 
#                 Spearman values for the same matrix location 
#                 (0.005 and 0.995 quantiles). The 1-1 line is also
#                 plotted and should pass between the quantiles for 
#                 the large majority of location pairs.
#
#Output
#A N by n by numsurrog array, surrogs. The ith surrogate is stored
#surrog[,,i].
#
ncsurrog<-function(m,numsurrog,plotcheckloc=NA)
{
  N<-dim(m)[1]
  n<-dim(m)[2]
  
  #Some error checking
  if (!all(is.finite(m)))
  {
    stop("Error in ncsurrog: NA, Inf, NaN, etc. not allowed in m")
  }
  
  #get spearman correlation matrix of the data
  scor<-cor(m,method="spearman")
    
  #get the covariance matrix of an mv normal with all
  #marginal variances 1 and with spearman correlation matrix 
  #scor
  ncov<-copula::iRho(copula::normalCopula(0,2),scor)
  ev<-eigen(ncov)$values
  
  #more error checking
  if(!matrixcalc::is.positive.semi.definite(ncov)){
    #stop("Error in ncsurrog: ncov is not positive semidefinite")
    warning("ncov is not +ve semi-definite, replacing with nearest +ve semi-definite matrix",call.=T,immediate.=T)

    ncovn<-Matrix::nearPD(ncov,corr=TRUE, maxit=1000)
    ncovn<-as.matrix(ncovn$mat)
    evn<-eigen(ncovn)$values
    comp_ev<<-cbind(ev,evn)
    colnames(comp_ev)<-c("evals_of_old_mat","evals_of_new_mat")
    cat("----This is comparing eigen values of the old and replacement matrices------","\n",comp_ev,"\n")
    ncov<-ncovn  # Now, replacing ncov by ncovn
  }
  
  #generate a bunch of mv normals in the shape of the final
  #desired output. Each row is a draw from an mv normal with 
  #mean 0 and cov matrix ncov.
  surrogs<-array(mvtnorm::rmvnorm(N*numsurrog,mean=rep(0,n),sigma=ncov),
                 dim=c(N,numsurrog,n))
  surrogs<-aperm(surrogs,c(1,3,2)) #now it is T by n by numsurrog
  
  #now align the ranks, see the alignranks function docs for details
  sm<-apply(FUN=sort,X=m,MARGIN=2)
  surrogs<-alignranks(sm,surrogs)
  
  # For additional check with optional plot
  if(!is.na(plotcheckloc))
  {
    #compute the Spearman correlation matrix for each surrogate
    scor_sur<-array(NA,c(n,n,numsurrog))
    for (counter in 1:numsurrog)
    {
      scor_sur[,,counter]<-cor(surrogs[,,counter],method="spearman")
    }
    
    #get quantiles of the distribution for each surrogate
    scor_quant<-apply(FUN=quantile,X=scor_sur,MARGIN=c(1,2),prob=c(0.005,0.995))
    scor_quant<-aperm(scor_quant,c(2,3,1))
    
    #plot the real values against the distributions of surrogate values
    scor<-scor[lower.tri(scor)]
    scor_lq<-scor_quant[,,1]
    scor_lq<-scor_lq[lower.tri(scor_lq)]
    scor_uq<-scor_quant[,,2]
    scor_uq<-scor_uq[lower.tri(scor_uq)]
    
    xlimits<-range(scor)
    ylimits<-range(scor_lq,scor_uq)
    
    panwd<-4
    panht<-4
    xht<-1
    ywd<-1
    gap<-.1
    totwd<-ywd+panwd+gap
    totht=xht+panht+gap
    grDevices::pdf(file=paste0(plotcheckloc,".pdf"),width=totwd,height=totht)
    par(fig=c(ywd/totwd,
              (ywd+panwd)/totwd,
              (xht)/totht,
              (xht+panht)/totht),
        mai=c(0,0,0,0),mgp=c(3,0.75,0))
    plot(c(-1,1),c(-1,1),type="l",lty="dashed",
         xlim=xlimits,ylim=ylimits)
    for (counter in 1:length(scor))
    {
      lines(rep(scor[counter],2),c(scor_lq[counter],scor_uq[counter]),type="l",col="grey")
    }
    mtext("Data correlation",side=1,line=2)
    mtext("Surrogate correlation quantiles",side=2,line=2)
    grDevices::dev.off()
  }

  return(surrogs)
}