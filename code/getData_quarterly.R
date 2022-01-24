## Download kelp data package from EDI and format it
## -- This is for annualized data.

rm(list=ls())


## 1. Download relevant data files ----------------------------------------------------------------

## kelp, nitrates, and waves
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/144/1/42e5ca2a0a1b51407aefb1a22cb61866" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")


dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "site_id",     
                 "year",     
                 "quarter",     
                 "kelp",     
                 "no3",     
                 "waves"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt1$site_id)!="factor") dt1$site_id<- as.factor(dt1$site_id)
if (class(dt1$quarter)!="factor") dt1$quarter<- as.factor(dt1$quarter)
if (class(dt1$kelp)=="factor") dt1$kelp <-as.numeric(levels(dt1$kelp))[as.integer(dt1$kelp) ]               
if (class(dt1$kelp)=="character") dt1$kelp <-as.numeric(dt1$kelp)
if (class(dt1$no3)=="factor") dt1$no3 <-as.numeric(levels(dt1$no3))[as.integer(dt1$no3) ]               
if (class(dt1$no3)=="character") dt1$no3 <-as.numeric(dt1$no3)
if (class(dt1$waves)=="factor") dt1$waves <-as.numeric(levels(dt1$waves))[as.integer(dt1$waves) ]               
if (class(dt1$waves)=="character") dt1$waves <-as.numeric(dt1$waves)


#climate indices
inUrl2  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/144/1/cfc9dd06fab598c27a78709c915ce81d" 
infile2 <- tempfile()
try(download.file(inUrl2,infile2,method="curl"))
if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")


climate_inds <-read.csv(infile2,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "year",     
                 "quarter",     
                 "NPGO",     
                 "MEI",     
                 "PDO"    ), check.names=TRUE)

unlink(infile2)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(climate_inds$quarter)!="factor") climate_inds$quarter<- as.factor(climate_inds$quarter)
if (class(climate_inds$NPGO)=="factor") climate_inds$NPGO <-as.numeric(levels(climate_inds$NPGO))[as.integer(climate_inds$NPGO) ]               
if (class(climate_inds$NPGO)=="character") climate_inds$NPGO <-as.numeric(climate_inds$NPGO)
if (class(climate_inds$MEI)=="factor") climate_inds$MEI <-as.numeric(levels(climate_inds$MEI))[as.integer(climate_inds$MEI) ]               
if (class(climate_inds$MEI)=="character") climate_inds$MEI <-as.numeric(climate_inds$MEI)
if (class(climate_inds$PDO)=="factor") climate_inds$PDO <-as.numeric(levels(climate_inds$PDO))[as.integer(climate_inds$PDO) ]               
if (class(climate_inds$PDO)=="character") climate_inds$PDO <-as.numeric(climate_inds$PDO)


## Coastline segment coordinates
inUrl5  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/144/1/37f8a10e2448275459b7055ca6e423cd" 
infile5 <- tempfile()
try(download.file(inUrl5,infile5,method="curl"))
if (is.na(file.size(infile5))) download.file(inUrl5,infile5,method="auto")


coords <-read.csv(infile5,header=F 
                  ,skip=1
                  ,sep=","  
                  ,quot='"' 
                  , col.names=c(
                    "site_id",     
                    "lat",     
                    "lon"    ), check.names=TRUE)

unlink(infile5)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(coords$site_id)!="factor") coords$site_id<- as.factor(coords$site_id)
if (class(coords$lat)=="factor") coords$lat <-as.numeric(levels(coords$lat))[as.integer(coords$lat) ]               
if (class(coords$lat)=="character") coords$lat <-as.numeric(coords$lat)
if (class(coords$lon)=="factor") coords$lon <-as.numeric(levels(coords$lon))[as.integer(coords$lon) ]               
if (class(coords$lon)=="character") coords$lon <-as.numeric(coords$lon)


## 2. Convert data matrices from long to wide format (suitable for analysis using the 'wsyn' package)

nsites <- length(unique(dt1$site_id))

dt1$yearquarter <- paste0(dt1$year, dt1$quarter)
ntimesteps <- length(unique(dt1$yearquarter))

site_ids <- unique(dt1$site_id)
yearquarters <- unique(dt1$yearquarter)

kelp <- matrix(NA, nsites, ntimesteps) #sites-by-years data matrix
no3 <- matrix(NA, nsites, ntimesteps)
waves <- matrix(NA, nsites, ntimesteps)

for(ii in 1:nsites){
  for(jj in 1:ntimesteps){
    kelp[ii,jj] <- dt1$kelp[dt1$site_id==site_ids[ii] & dt1$yearquarter==yearquarters[jj]]
    no3[ii,jj] <- dt1$no3[dt1$site_id==site_ids[ii] & dt1$yearquarter==yearquarters[jj]]
    waves[ii,jj] <- dt1$waves[dt1$site_id==site_ids[ii] & dt1$yearquarter==yearquarters[jj]]
  }
}
