## function for plotting segments perpendicular to the coastline

coastsegs<-function(x0, y0, coast.coords, len, lwd=1, color){
  
  if(length(x0)!=length(y0)){
    stop("x0 and y0 must be the same length")
  }
  
  nn <- length(x0)
  
  for(ii in 1:nn){
    
    #find nearest pair of coastline points 
  
    #calculate distances
    dd <- NULL
    for(jj in 1:nrow(coast.coords)){
      
      dd <- c(dd, sqrt( (x0[ii]-coast.coords$X[jj])^2 + (y0[ii]-coast.coords$Y[jj])^2))
      
    }
    
    xn <- coast.coords$X[which.min(dd)]
    yn <- coast.coords$Y[which.min(dd)]
    
    a <- atan2(yn-y0[ii], xn-x0[ii])
    
    x1 = x0[ii] + len * cos(a + pi/2)
    y1 = y0[ii] + len * sin(a + pi/2)
    
    segments(x0[ii], y0[ii], x1, y1, lwd=lwd, col=color[ii])
  }
  
  
}