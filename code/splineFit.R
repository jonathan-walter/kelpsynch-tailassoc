## Spline correlograms for matrices of pairwise wavelet coherence!
##  Modification of the 'spline.correlog' fucntion from package ncf and for full functionality 'ncf' should be active in the R session
##  Note: this will also work on any user-input matrix of similarities between locations, regardless of what method produced it
## Inputs:
##  coords = n by 2 matrix, where n is the number of sites and columns are x and y coordinates
##  sim.mat = an n by n matrix of similarities 
##  df =  spline degrees of freedom, defaults to sqrt(n)
##  boot.resamp = number of resamples for bootstrapped confidence intervals
##  npoints = the number of points at which to save the value for the spline function (and confidence intervals)
##  save = if TRUE the whole matrix of output from resampling is saved as an resamp by npoints dimensional matrix
##  filter = if TRUE the Fourier filter method of Hall et al. is applied to ensure positive semidifinitenss of the estimator
##  fw = if filter is TRUE, fw may be used to truncate the function at some distance given by fw. If fw = 0 (the default), no truncation is done.
##  max.it = maximum number of iterations for calculating the intercepts using the Newton method
##  xmax = maximum distance used. If FALSE, the maximum in the observed data is used.
##  latlon = if TRUE, coordinates are latitude and longitude, and 'gcdist' in package 'ncf' is used to calculate the spatial distance (km) from
##    lat-long data. This is not recommended as others have found this to give unexpected results. Rather, it is better to first project
##    coordinates into a grid-based coordinate system. If distances are small enough (e.g., <200 km) it is probably OK to treat lat-long 
##    coordinates as a grid.
##  quiet =  if TRUE the resampling counter is suppressed during execution.


###### For testing during development, uncomment to use or delete
# set.seed(11)
# n = 20
# sim.mat<-matrix(runif(n*n,0,1), nrow=n, ncol=n)
# coords<-cbind(runif(n,0,50),runif(n,0,50))
# #df = NULL; boot.resamp = 10; npoints = 300; save = FALSE; filter = FALSE; fw = 0; max.it = 25; xmax = FALSE; latlon = FALSE; quiet = TRUE
# 
# # test<-spline.correlog.simmat(coords, sim.mat)
# # plot(test)

## start the function!
splineFit<-function(coords, sim.mat, df=NULL, boot.resamp=1000, npoints=300, save=FALSE, 
                              filter=FALSE, fw=0, max.it=25, xmax=FALSE, latlon=FALSE, quiet=TRUE)
{
  real <- list(x.intercept = NA, e.intercept = NA, y.intercept = NA, 
               predicted = list(x = matrix(NA, nrow = 1, ncol = npoints), 
                                y = matrix(NA, nrow = 1, ncol = npoints))) ##prepare output
  n=nrow(sim.mat)
  if (is.null(df)){df <- sqrt(n)}
  x = coords[,1]
  y = coords[,2]
  ## Calculate 'xdist' -- the distance metric
  if (latlon) {
    xdist <- gcdist(x,y)
  }
  else {
    xdist <- sqrt(outer(x, x, "-")^2 + outer(y, y, "-")^2) #distances based on x and y coordinates
  }
  maxdist <- ifelse(!xmax, max(na.omit(xdist)), xmax)
  ## Get ready to compute empirical ("real") spline
  triang <- lower.tri(xdist)
  u <- xdist[triang] #x in smooth.spline
  v <- sim.mat[triang] #y in smooth.spline -- this is our relatedness metric
  sel <- is.finite(u) & is.finite(v)
  u <- u[sel]
  v <- v[sel]
  v <- v[u <= maxdist]
  u <- u[u <= maxdist]
  sobj <- smooth.spline(u, v, df = df) # this is the empirical smoothing spline
  xpoints <- seq(0, maxdist, length = npoints)
  lx <- predict(sobj, x = xpoints)
  if (filter == TRUE) {
    if (fw > 0) {
      lx$y[xpoints > fw] <- 0
    }
    lx$y <- ff.filter(lx$y)
  }
  real$y.intercept <- lx$y[1]
  real$predicted <- list(x = xpoints, y = lx$y)
  konst <- 1  # Calculating the x-intercept distance
  if (real$y.intercept < 0) { #switch things around if spline starts negative
    lx$y <- -lx$y 
    konst <- -1
  }
  ly <- 1:length(lx$y)
  choise <- ly[lx$y < 0][1] #find first negative spline y-value
  pos <- lx$x[choise - 1]
  neg <- lx$x[choise]
  pos <- pos + (neg - pos)/2
  tmp <- smooth.spline(lx) 
  for (j in 1:max.it) { #this is iterating through something, I believe to refine the value of x.intercept
    if (is.na(neg)) {
      pos <- NA
      break
    }
    if (neg == 0) {
      pos <- 0
      break
    }
    neg <- pos - predict(tmp, pos)$y/predict(tmp, pos, deriv = 1)$y
    if (abs(pos - neg) < 1e-06) {
      break
    }
    pos <- neg
  }
  real$x.intercept <- konst * pos #done calculating the value of the x.intercept
  sobj <- smooth.spline(u, v - 1/exp(1), df = df) # calculating the value of e.intercept
  lx <- predict(sobj, x = xpoints)                # which I think has meaning if your data
  if (filter == TRUE) {                           # are ln(y + 1) transformed
    if (fw > 0) {
      lx$y[xpoints > fw] <- 0
    }
    lx$y <- ff.filter(lx$y)
  }
  ly <- 1:length(lx$y)
  choise <- ly[lx$y < 0][1]
  pos <- lx$x[choise - 1]
  neg <- lx$x[choise]
  pos <- pos + (neg - pos)/2
  tmp <- smooth.spline(lx)
  for (j in 1:max.it) {
    if (is.na(neg)) {
      pos <- NA
      break
    }
    if (neg == 0) {
      pos <- 0
      break
    }
    neg <- pos - predict(tmp, pos)$y/predict(tmp, pos, deriv = 1)$y
    if (abs(pos - neg) < 1e-06) {
      break
    }
    pos <- neg
  }
  real$e.intercept <- pos
  boot <- list(NULL) ## Start bootstrapping!
  boot$boot.summary <- list(NULL)
  if (boot.resamp != 0) {
    boot$boot.summary$x.intercept <- matrix(NA, nrow = boot.resamp, ncol = 1)
    boot$boot.summary$e.intercept <- matrix(NA, nrow = boot.resamp, ncol = 1)
    boot$boot.summary$y.intercept <- matrix(NA, nrow = boot.resamp, ncol = 1)
    predicted <- list(x = matrix(NA, nrow = 1, ncol = npoints), y = matrix(NA, nrow = boot.resamp, ncol = npoints))
    predicted$x[1, ] <- xpoints
#     type <- charmatch(type, c("boot", "perm"), nomatch = NA)
#     if (is.na(type)) 
#       stop("method should be \"boot\", or \"perm\"")
    for (i in 1:boot.resamp) {
      if (!quiet) {cat(i, " of ", boot.resamp, "\n")}
      n=nrow(sim.mat)
      trekkx <- sample(1:n, replace = TRUE)
      trekky <- trekkx
      
      xdistb <- xdist[trekkx, trekkx] ## Permuting the locations (distance part)
      triang <- lower.tri(xdist)
      xdistb <- xdistb[triang]
      sim.matb <- sim.mat[trekky, trekky][triang] ## Permuting the locations (correlation part)

      sim.matb <- sim.matb[!(xdistb == 0)]
      xdistb <- xdistb[!(xdistb == 0)]

      u <- xdistb
      v <- sim.matb
      sel <- is.finite(u) & is.finite(v)
      u <- u[sel]
      v <- v[sel]
      v <- v[u <= maxdist]
      u <- u[u <= maxdist]
      sobj <- smooth.spline(u, v, df = df)
      lx <- predict(sobj, x = xpoints)
      if (filter == TRUE) {
        if (fw > 0) {
          lx$y[xpoints > fw] <- 0
        }
        lx$y <- ff.filter(lx$y)
      }
      boot$boot.summary$y.intercept[i, 1] <- lx$y[1]
      predicted$y[i, ] <- lx$y
      konst <- 1
      if (boot$boot.summary$y.intercept[i, 1] < 0) {
        lx$y <- -lx$y
        konst <- -1
      }
      ly <- 1:length(lx$y)
      choise <- ly[lx$y < 0][1]
      pos <- lx$x[choise - 1]
      neg <- lx$x[choise]
      pos <- pos + (neg - pos)/2
      tmp <- smooth.spline(lx)
      for (j in 1:max.it) {
        if (is.na(neg)) {
          pos <- NA
          break
        }
        if (neg == 0) {
          pos <- 0
          break
        }
        neg <- pos - predict(tmp, pos)$y/predict(tmp, 
                                                 pos, deriv = 1)$y
        if (abs(pos - neg) < 1e-06) {
          break
        }
        pos <- neg
      }
      boot$boot.summary$x.intercept[i, 1] <- konst * pos
      sobj <- smooth.spline(u, v - 1/exp(1), df = df)
      lx <- predict(sobj, x = xpoints)
      if (filter == TRUE) {
        if (fw > 0) {
          lx$y[xpoints > fw] <- 0
        }
        lx$y <- ff.filter(lx$y)
      }
      ly <- 1:length(lx$y)
      choise <- ly[lx$y < 0][1]
      pos <- lx$x[choise - 1]
      neg <- lx$x[choise]
      pos <- pos + (neg - pos)/2
      tmp <- smooth.spline(lx)
      for (j in 1:max.it) {
        if (is.na(neg)) {
          pos <- NA
          break
        }
        if (neg == 0) {
          pos <- 0
          break
        }
        neg <- pos - predict(tmp, pos)$y/predict(tmp, 
                                                 pos, deriv = 1)$y
        if (abs(pos - neg) < 1e-06) 
          break
        pos <- neg
      }
      boot$boot.summary$e.intercept[i, 1] <- pos
    }
    if (save == TRUE) {
      boot$boot <- list(predicted = predicted)
    }
    else {
      boot$boot <- NULL
    }
    ty <- apply(predicted$y, 2, quantile, probs = c(0, 0.025, 
                                                    0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1), 
                na.rm = TRUE)
    dimnames(ty) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 
                           0.75, 0.9, 0.95, 0.975, 1), NULL)
    tx <- predicted$x
    boot$boot.summary$predicted <- list(x = tx, y = ty)
  }
  else {
    boot <- NULL
    boot.summary <- NULL
  }
  res <- list(real = real, boot = boot, max.distance = maxdist, 
              call = deparse(match.call()))
  class(res) <- "spline.correlog"
  return(res)
}