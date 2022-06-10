# Adapt Colias niche model

#Center for urban horticulture
#lat 47.657628, Longitude: -122.290255

#Corfu site in central Washington in the Columbia National Wildlife Refuge
#12 miles west of Othello
#46.850663264 -119.535331192

#------------------
#setwd("/Users/laurenbuckley/Downloads/")

library(NicheMapR)
library(ecmwfr)
library(mcera5)
library(lubridate)
library(dplyr)
library(tidync)

#install.packages("remotes")
#remotes::install_github("everydayduffy/mcera5")

# get ERA5 data with package mcera5 (just do once for region and time of interest)
# assign your credentials (register here: https://cds.climate.copernicus.eu/user/register)
uid <- "78176"
cds_api_key <- "062c3e77-bcc8-4c56-8e72-4872e7a92be6"

ecmwfr::wf_set_key(user = uid, key = cds_api_key, service = "cds")

#locations Corfu, P. occidentalis, Apr to Sept
#Seattle, P. rapae,
locations= c("Corfu","Seattle")

for(loc.k in 1:2){

# bounding coordinates (in WGS84 / EPSG:4326)
if(loc.k==1){xmn <- -119.6; xmx <- -119.45; ymn <- 46.75; ymx <- 47}
if(loc.k==2){xmn <- -122.45; xmx <- -122.2; ymn <- 47.5; ymx <- 47.75}

if(loc.k==1){loc <- c(-119.535331192, 46.850663264)}  #Corfu
if(loc.k==2){loc <- c(-122.290255, 47.657628)} #Seattle  
  
#years for data
#if(loc.k==1) years=c(1989:1993, 2017:2021)
#if(loc.k==2) years=c(2001:2005, 2017:2021)
if(loc.k==1) years=c(2002:2017)
if(loc.k==2) years=c(2009:2016) ##Error with 2008

#set microclim path
file_prefix="era5"
if(loc.k==1) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/',file_prefix, sep="")
if(loc.k==2) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/',file_prefix, sep="")
# filename and location for downloaded .nc files
if(loc.k==1) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/'
if(loc.k==2) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/'

#check works

for(yr in years){

#set dates for era5 call and NicheMapr
st_date<- paste(yr,":04:01",sep="")  
en_date<- paste(yr,":09:30",sep="") 

dstart <- paste("01/04/", yr,sep="")
dfinish <- paste("30/09/", yr,sep="")

# temporal extent
st_time <- lubridate::ymd(st_date)
en_time <- lubridate::ymd(en_date)

# build a request (covering multiple years)
req <- build_era5_request(xmin = xmn, xmax = xmx,
                          ymin = ymn, ymax = ymx,
                          start_time = st_time,
                          end_time = en_time,
                          outfile_name = file_prefix)
str(req)
request_era5(request = req, uid = uid, out_path = op)

# run micro_era5 for a location (make sure it's within the bounds of your .nc files)
micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=0.01, runshade = 1, spatial = spatial_path)
#https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html
#minshade, Minimum shade level to use (can be a single value or a vector of daily values) (%)
#maxshade, Maximum shade level to use (can be a single value or a vector of daily values) (%)
#Usrhyt, Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
# runshade = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)

#test with direct download
#file_prefix="Corfu"
#spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5/',file_prefix, sep="")

#--------------

#combine and save data 
metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade

# append dates
tzone<-paste("Etc/GMT+",0,sep="")
dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")

metout <- cbind(dates,metout)
soil <- cbind(dates,soil)
#soilmoist <- cbind(dates, soilmoist)

#extract air temperature and soil
#dat <- cbind(metout[,1:3],metout$TAREF, metout$ZEN, metout$SOLR, soil$D0cm)
dat <- cbind(metout, soil[,4:13])

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro/")
write.csv(dat, paste(locations[loc.k], yr,".csv",sep=""))

} #end loop years
} #end loop locations

#-----------------------
#read back in
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro/")
metout= read.csv("Seattle2020.csv")
#fix dates
metout$dates= as.Date(metout$dates)

# plotting above-ground conditions in minimum shade
#temp
with(metout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                  , type = "l",main=paste("air and sky temperature",sep=""), ylim = c(-20, 60))})
with(metout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                    , type = "l",lty=2,col='blue')})
with(metout,{points(TSKYC ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                    ,  type = "l",col='light blue',main=paste("sky temperature",sep=""))})

with(metout,{plot(RHLOC ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
                  , type = "l",ylim=c(0,100),main=paste("humidity",sep=""))})
with(metout,{points(RH ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
                    , type = "l",col='blue',lty=2,ylim=c(0,100))})
with(metout,{plot(VREF ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
                  ,  type = "l",main="wind speed",ylim = c(0, 15))})
with(metout,{points(VLOC ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
                    ,  type = "l",lty=2,col='blue')})
with(metout,{plot(SOLR ~ dates,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
                  ,  type = "l",main="solar radiation")})
with(metout,{plot(SNOWDEP ~ dates,xlab = "Date and Time", ylab = "Snow Depth (cm)"
                  ,  type = "l",main="snow depth")})

# plotting soil temperature
for(i in 1:10){
  if(i==1){
    plot(soil[,i+3]~soil[,1],xlab = "Date and Time", ylab = "Soil Temperature (Â°C)"
         ,col=i,type = "l",main=paste("soil temperature",sep=""))
  }else{
    points(soil[,i+3]~soil[,1],xlab = "Date and Time", ylab = "Soil Temperature
    (Â°C)",col=i,type = "l")
  }
}

# plotting soil moisture
for(i in 1:10){
  if(i==1){
    plot(soilmoist[,i+3]*100~soilmoist[,1],xlab = "Date and Time", ylab = "Soil Moisture (% volumetric)"
         ,col=i,type = "l",main=paste("soil moisture",sep=""))
  }else{
    points(soilmoist[,i+3]*100~soilmoist[,1],xlab = "Date and Time", ylab = "Soil Moisture
    (%)",col=i,type = "l")
  }
}

#-----------------------
#Check microclimate analysis

#locations Corfu, P. occidentalis, Apr to Sept
#Seattle, P. rapae,
locations= c("Corfu","Seattle")

loc.k<- 1

if(loc.k==1){loc <- c(-119.535331192, 46.850663264)}  #Corfu
if(loc.k==2){loc <- c(-122.290255, 47.657628)} #Seattle  

  if(loc.k==1) years=c(1989:2002,2018:2021) #1989:2021
  if(loc.k==2) years=c(2009:2016) ##Error with 2008
  
  #set microclim path
  file_prefix="era5"
  if(loc.k==1) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/',file_prefix, sep="")
  if(loc.k==2) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/',file_prefix, sep="")
  # filename and location for downloaded .nc files
  if(loc.k==1) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/'
  if(loc.k==2) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/'
  
  for(yr in years){
    
    #set dates for era5 call and NicheMapr
    dstart <- paste("01/04/", yr,sep="")
    dfinish <- paste("30/09/", yr,sep="")
    
    # run micro_era5 for a location (make sure it's within the bounds of your .nc files)
    micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=0.5, runshade = 0, spatial = spatial_path, minshade=0, maxshade=10)
    #micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=0.5, runshade = 0, spatial = spatial_path, minshade=90, maxshade=100)
    
    #https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html
    #minshade, Minimum shade level to use (can be a single value or a vector of daily values) (%)
    #maxshade, Maximum shade level to use (can be a single value or a vector of daily values) (%)
    #Usrhyt, Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
    # runshade = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)
    
    #--------------
    
    #combine and save data 
    metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
    soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
    soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
    
    # append dates
    tzone<-paste("Etc/GMT+",0,sep="")
    dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
    
    metout <- cbind(dates,metout)
    soil <- cbind(dates,soil)
    #soilmoist <- cbind(dates, soilmoist)
    
    #extract air temperature and soil
    #dat <- cbind(metout[,1:3],metout$TAREF, metout$ZEN, metout$SOLR, soil$D0cm)
    dat <- cbind(metout, soil[,4:13])
    
    setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/")
    #setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_shade/")
    
    write.csv(dat, paste(locations[loc.k], yr,".csv",sep=""))
    
  } #end loop years


