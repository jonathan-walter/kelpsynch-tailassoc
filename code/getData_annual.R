## Download kelp data package from EDI and format it
## -- This is for annualized data.

rm(list=ls())


## 1. Download relevant data files ----------------------------------------------------------------

## kelp biomass, nitrate concentration, and wave height
inUrl3  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/144/1/2aa452e37d7a6ce03e8af19311d924dd" 
infile3 <- tempfile()
try(download.file(inUrl3,infile3,method="curl"))
if (is.na(file.size(infile3))) download.file(inUrl3,infile3,method="auto")


dt3 <-read.csv(infile3,header=F 
               ,skip=1
               ,sep=","  
               , col.names=c(
                 "site_id",     
                 "year",     
                 "kelp",     
                 "no3",     
                 "waves"    ), check.names=TRUE)

unlink(infile3)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

if (class(dt3$site_id)!="factor") dt3$site_id<- as.factor(dt3$site_id)
if (class(dt3$kelp)=="factor") dt3$kelp <-as.numeric(levels(dt3$kelp))[as.integer(dt3$kelp) ]               
if (class(dt3$kelp)=="character") dt3$kelp <-as.numeric(dt3$kelp)
if (class(dt3$no3)=="factor") dt3$no3 <-as.numeric(levels(dt3$no3))[as.integer(dt3$no3) ]               
if (class(dt3$no3)=="character") dt3$no3 <-as.numeric(dt3$no3)
if (class(dt3$waves)=="factor") dt3$waves <-as.numeric(levels(dt3$waves))[as.integer(dt3$waves) ]               
if (class(dt3$waves)=="character") dt3$waves <-as.numeric(dt3$waves)


## Climate indices
inUrl4  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-sbc/144/1/91190eac2b6b35f58353a19dbba676b9" 
infile4 <- tempfile()
try(download.file(inUrl4,infile4,method="curl"))
if (is.na(file.size(infile4))) download.file(inUrl4,infile4,method="auto")


climate_inds <-read.csv(infile4,header=F 
               ,skip=1
               ,sep=","  
               ,quot='"' 
               , col.names=c(
                 "year",     
                 "NPGO",     
                 "MEI",     
                 "PDO"    ), check.names=TRUE)

unlink(infile4)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

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

nsites <- length(unique(dt3$site_id))
nyears <- length(unique(dt3$year))

site_ids <- unique(dt3$site_id)
years <- unique(dt3$year)

kelp <- matrix(NA, nsites, nyears) #sites-by-years data matrix
no3 <- matrix(NA, nsites, nyears)
waves <- matrix(NA, nsites, nyears)

for(ii in 1:nsites){
  for(jj in 1:nyears){
    kelp[ii,jj] <- dt3$kelp[dt3$site_id==site_ids[ii] & dt3$year==years[jj]]
    no3[ii,jj] <- dt3$no3[dt3$site_id==site_ids[ii] & dt3$year==years[jj]]
    waves[ii,jj] <- dt3$waves[dt3$site_id==site_ids[ii] & dt3$year==years[jj]]
  }
}
