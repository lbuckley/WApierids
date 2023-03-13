#Envi data to process
#Sacramento: shade 2011-2016, sun 1961:1971, 2001:2021
#Los Banos: run sun and shade

# Adapt Colias niche model

#Center for urban horticulture
#lat 47.657628, Longitude: -122.290255

#Corfu site in central Washington in the Columbia National Wildlife Refuge
#12 miles west of Othello
#46.850663264 -119.535331192

#Nielsen and Kingsolver Site: Isleton, California, United States (38.154°N, 121.677°W, 0 m)
#Initial Hoffman site: Los Banos, California (37.06°N, 120.85°W, 36 m)
#Initial Hoffman 1978: https://www.jstor.org/stable/pdf/2460345.pdf

#Higgins
#Sacramento Valley: 38.44N, 121.86W
#Montrose CO: 38.62N, 108.02W
#Initial Sherman and Watt: https://link.springer.com/article/10.1007/BF00694570
#Los Banos, Merced Co., California, and from Hotchkiss, Delta Co., Colorado

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
locations= c("Corfu","Seattle","Montrose","Sacramento","LosBanos")

#add additional TPCs
#Higgins: 
#Montrose Valley, CO: 1961-1971; 2001-2011;N 38.62, W 108.02, 1633m
#Sacramento Valley, CA: 1961-1971; 2001-2011; N 38.44, W121.86, 19m

#Nielsen
#Los Banos, CA: ; 37.06N, 120.85W, 36M, iniital 1970/1971, 2018;
#data: 1962-1971, 2009-2018

for(loc.k in 6){ #3:5

# bounding coordinates (in WGS84 / EPSG:4326)
if(loc.k==1){xmn <- -119.6; xmx <- -119.45; ymn <- 46.75; ymx <- 47}
if(loc.k==2){xmn <- -122.45; xmx <- -122.2; ymn <- 47.5; ymx <- 47.75}
if(loc.k==3){xmn <- -108.03; xmx <- -108.01; ymn <- 38.61; ymx <- 38.63}
if(loc.k==4){xmn <- -121.87; xmx <- -121.85; ymn <- 38.43; ymx <- 38.45}
if(loc.k==5){xmn <- -120.86; xmx <- -120.84; ymn <- 37.05; ymx <- 37.07}
if(loc.k==6){xmn <- -121.87; xmx <- -121.85; ymn <- 38.43; ymx <- 38.45} #Sacramento for Nielsen pupation estimate
   
if(loc.k==1){loc <- c(-119.535331192, 46.850663264)}  #Corfu
if(loc.k==2){loc <- c(-122.290255, 47.657628)} #Seattle  
if(loc.k==3){loc <- c(-108.02,38.62)} #Montrose  
if(loc.k==4){loc <- c(-121.86, 38.44)} #Sacramento 
if(loc.k==5){loc <- c(-120.85, 37.06)} #LosBanos 
if(loc.k==6){loc <- c(-121.86, 38.44)} #Sacramento 
  
#years for data
#if(loc.k==1) years=c(1989:1993, 2017:2021)
#if(loc.k==2) years=c(2001:2005, 2017:2021)
if(loc.k==1) years=c(2002:2017)
if(loc.k==2) years=c(1995:2000) ##Error with 2008
if(loc.k==3) years= c(1972:2001)   #c(1961:1971, 2001:2011, 2012:2021)
if(loc.k==4) years= c(2017:2021)  #c(2017:2021,1962)  #c(1961:1971, 2001:2011, 2012:2021)
if(loc.k==5) years= c(1961:2021)
if(loc.k==6) years= c(2018) #c(1970,2018)

#set microclim path
file_prefix="era5"
if(loc.k==1) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/',file_prefix, sep="")
if(loc.k==2) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/',file_prefix, sep="")
if(loc.k==3) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Montrose/',file_prefix, sep="")
if(loc.k==4) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento/',file_prefix, sep="")
if(loc.k==5) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_LosBanos/',file_prefix, sep="")
if(loc.k==6) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento_Nielsen/',file_prefix, sep="")

# filename and location for downloaded .nc files
if(loc.k==1) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/'
if(loc.k==2) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/'
if(loc.k==3) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Montrose/'
if(loc.k==4) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento/'
if(loc.k==5) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_LosBanos/'
if(loc.k==6) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento_Nielsen/'

for(yr in years){ #years

#set dates for era5 call and NicheMapr
st_date<- paste(yr,":04:01",sep="") 
if(loc.k %in% c(4,5)) st_date<- paste(yr,":03:01",sep="")
if(loc.k %in% c(6)) st_date<- paste(yr,":01:01",sep="")
en_date<- paste(yr,":09:30",sep="") 

dstart <- paste("01/04/", yr,sep="")
if(loc.k %in% c(4,5)) dstart <- paste("01/03/", yr,sep="")
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

# # run micro_era5 for a location (make sure it's within the bounds of your .nc files)
# micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=0.2, runshade = 0, spatial = spatial_path)
# #https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html
# #minshade, Minimum shade level to use (can be a single value or a vector of daily values) (%)
# #maxshade, Maximum shade level to use (can be a single value or a vector of daily values) (%)
# #Usrhyt, Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
# # runshade = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)
# 
# #test with direct download
# #file_prefix="Corfu"
# #spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5/',file_prefix, sep="")
# 
# #--------------
# 
# #combine and save data 
# metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
# soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
# soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
# 
# # append dates
# tzone<-paste("Etc/GMT+",0,sep="")
# dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
# 
# metout <- cbind(dates,metout)
# soil <- cbind(dates,soil)
# #soilmoist <- cbind(dates, soilmoist)
# 
# #extract air temperature and soil
# #dat <- cbind(metout[,1:3],metout$TAREF, metout$ZEN, metout$SOLR, soil$D0cm)
# dat <- cbind(metout, soil[,4:13])
# 
# setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro/")
# write.csv(dat, paste(locations[loc.k], yr,".csv",sep=""))

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

#============================
#Check microclimate analysis

#locations Corfu, P. occidentalis, Apr to Sept
#Seattle, P. rapae,
locations= c("Corfu","Seattle","Montrose","Sacramento","LosBanos")

loc.k<- 2 #2,3,4
shade<- F

height= c(0.3,0.2,0.2,0.2,0.2,0.2)[loc.k]

if(loc.k==1){loc <- c(-119.535331192, 46.850663264)}  #Corfu
if(loc.k==2){loc <- c(-122.290255, 47.657628)} #Seattle  
if(loc.k==3){loc <- c(-108.02,38.62)} #Montrose  
if(loc.k==4){loc <- c(-121.86, 38.44)} #Sacramento 
if(loc.k==5){loc <- c(-120.85, 37.06)} #Los Banos
if(loc.k==6){loc <- c(-121.86, 38.44)} #Sacramento 

  if(loc.k==1) years=c(1989:2021) 
  if(loc.k==2) years=c(1998:2021) ##2001:2021
  if(loc.k==3) years=c(1961:2011)  #c(1961:1971, 2001:2011) 
  if(loc.k==4) years=c(2018:2021)  #c(1971,2001:2021) #c(1961:1968,1970:1971,2001:2021) #c(1961:1968,1970:1971,2001:2021)  #c(1961:1971, 2001:2016)
  if(loc.k==5) years=c(1962:1971, 2009:2018)  #c(1961:1971, 2001:2016)
  if(loc.k==6) years= c(1970,2018)

  #set microclim path
  file_prefix="era5"
  if(loc.k==1) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/',file_prefix, sep="")
  if(loc.k==2) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/',file_prefix, sep="")
  if(loc.k==3) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Montrose/',file_prefix, sep="")
  if(loc.k==4) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento/',file_prefix, sep="")
  if(loc.k==5) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_LosBanos/',file_prefix, sep="")
  if(loc.k==6) spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento_Nielsen/',file_prefix, sep="")
  
  
  # filename and location for downloaded .nc files
  if(loc.k==1) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Corfu/'
  if(loc.k==2) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Seattle/'
  if(loc.k==3) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Montrose/'
  if(loc.k==4) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento/'
  if(loc.k==5) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_LosBanos/'
  if(loc.k==6) op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento_Nielsen/'
  
  for(yr in years){
    
    #set dates for era5 call and NicheMapr
    dstart <- paste("01/04/", yr,sep="")
    if(loc.k==6) dstart <- paste("01/01/", yr,sep="")
    dfinish <- paste("30/09/", yr,sep="")
    
    # run micro_era5 for a location (make sure it's within the bounds of your .nc files)
    if(shade==F) micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=height, runshade = 0, spatial = spatial_path, minshade=25, maxshade=30, REFL=0.3, RUF=0.02, IR=2)
    if(shade==T) micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=height, runshade = 0, spatial = spatial_path, minshade=90, maxshade=100, REFL=0.3, RUF=0.02, IR=2)
    
    #https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html
    #https://github.com/mrke/NicheMapR/blob/master/R/micro_era5.R
    #minshade, Minimum shade level to use (can be a single value or a vector of daily values) (%)
    #maxshade, Maximum shade level to use (can be a single value or a vector of daily values) (%)
    #Usrhyt, Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
    # runshade = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)
    
    #' @param REFL Soil solar reflectance, decimal \%
    #' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1) or from ERA5 data (2)}\cr\cr
    #' \code{RUF}{ = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.}\cr\cr
    
    #REFL= 0.3 #30% substrate solar reflectivity, Kingsolver 1983
    #RUF=0.02 is surface roughness height
    #--------------
    
    #combine and save data 
    metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
    soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
    soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
    #shadmet The above ground micrometeorological conditions under the maximum specified shade
    
    # append dates
    tzone<-paste("Etc/GMT+",0,sep="")
    dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
    
    metout <- cbind(dates,metout)
    soil <- cbind(dates,soil)
    #soilmoist <- cbind(dates, soilmoist)
    
    #extract air temperature and soil
    #dat <- cbind(metout[,1:3],metout$TAREF, metout$ZEN, metout$SOLR, soil$D0cm)
    dat <- cbind(metout, soil[,4:13])
    
    if(shade==F) setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/")
    if(shade==T) setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_shade/")
    
    write.csv(dat, paste(locations[loc.k], yr,".csv",sep=""))
   # write.csv(dat, paste(locations[loc.k], yr,"_Jan.csv",sep=""))
  } #end loop years

#library(raster)
#r= raster('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Sacramento/era5_1961.nc')
  
#======================
#plot output
  
  #plot temperature distributions
  locations= c("Montrose","Sacramento", "LosBanos")
  
 loc.k=1
    
    #years for data
 if (loc.k==1) years=c(1972:2014, 2016:2021) 
    #check 2015
    if (loc.k==2) years=c(1961:2021) 
    if (loc.k==3) years=c(1962:1971, 2009:2018)
    
    ##SUN
    setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/')
    #SHADE
    #setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_shade/')
    
    #combine data
    for(yr.k in 1:length(years)){
      dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
      
      
      dat$year= years[yr.k]
      
      dat$period<- NA
      #vary periods
      #if(years[yr.k]<2001)dat$period="initial"
      #if(years[yr.k]>=2001 & years[yr.k]<2011) dat$period="middle"
      #if(years[yr.k]>=2011)dat$period="recent"
      if(years[yr.k]<1972) dat$period<- "1961-1971"
      if(years[yr.k]>=2001 & years[yr.k]<=2011) dat$period <- "2001-2011"
      
      if(yr.k==1) dat.all=dat
      if(yr.k>1) dat.all=rbind(dat, dat.all)
    }
    
    #keep all hours, no longer subset to sunlight
    dat.day=dat.all #subset(dat.all, dat.all$SOLR>0)
    
    #Recode as Apr + May, July +Aug; 91:151, 182:243
    dat.day$seas= NA
    #April to May
    dat.day$seas[dat.day$DOY %in% 91:151]= "AprMay"
    # July to August
    dat.day$seas[dat.day$DOY %in% 182:243]= "JulAug"
    # #order
    dat.day$seas= factor(dat.day$seas, levels=c("AprMay","JulAug") )
    #drop other days
    dat.day= dat.day[which(!is.na(dat.day$seas)),]
    
    #combine doy and time
    dat.day$d.hr= dat.day$DOY + dat.day$TIME/60
    
    #drop middle period
    dat.day.plot= dat.day[!is.na(dat.day$period),]
    #make character factor
    dat.day.plot$period= as.factor(as.character(dat.day.plot$period))
    
    #--------------------------  
    #plot density distributions
    p1= ggplot(dat.day.plot, aes(x=TALOC))+
      geom_density(alpha=0.3, aes(fill=period, color=period))+
      facet_wrap(~seas)+
      ylab("Feeding rate (g/g/h)")+
      xlab("Temperature (°C)" )+
      xlim(-5,50)+
      theme_classic(base_size = 20) +
      theme(legend.position = "none")+
      scale_y_continuous(limits=c(0,0.102), expand=c(0,0))+
      scale_color_manual(values=c("#3CBB75FF","#404788FF")  )+
      scale_fill_manual(values=c("#3CBB75FF","#404788FF")  ) 
    #D0cm, TALOC, TAREF
    
    #plot reference temperatures
    p1.ref= ggplot(dat.day.plot, aes(x=TAREF))+
      geom_density(alpha=0.3, aes(fill=period, color=period))+
      facet_wrap(~seas)+
      ylab("Feeding rate (g/g/h)")+
      xlab("Temperature at reference height (°C)" )+
      xlim(-5,50)+ 
      theme_classic(base_size = 20) +
      theme(legend.position = c(0.4, 0.8),legend.background = element_rect(fill="transparent"))+
      scale_y_continuous(limits=c(0,0.102), expand=c(0,0))+
      scale_color_manual(values=c("#3CBB75FF","#404788FF")  )+
      scale_fill_manual(values=c("#3CBB75FF","#404788FF")  ) 
    
    #plot together
    p1.pl.ref=p1 + geom_density(alpha=0.3, linetype="dashed", aes(x=TAREF, color=period))+
      theme(legend.position = c(0.4, 0.8))
    
    #remove label
    p1= p1+theme(strip.text.x = element_blank())
    
    
    if(loc.k==1) {temp.co= p1; temp.co.ref= p1.ref; dat.day.co=dat.day; temps.co= p1.pl.ref}
    if(loc.k==2) {temp.ca= p1; temp.ca.ref= p1.ref; dat.day.ca=dat.day; temps.ca= p1.pl.ref; temp.nielsen=p1n.pl.ref}
    
#=========================
#check microclimate model
    
    #set microclim path
    file_prefix="era5"
   spatial_path<- paste('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Montrose/',file_prefix, sep="")
   op<- '/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_Montrose/'
   
   yr= 2010
   
      #set dates for era5 call and NicheMapr
      dstart <- paste("01/04/", yr,sep="")
      dfinish <- paste("30/09/", yr,sep="")
      
      par(mfrow=c(1,2))
      
      # run micro_era5 for a location (make sure it's within the bounds of your .nc files)
      micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=height, runshade = 0, spatial = spatial_path, minshade=0, maxshade=10, REFL=0.3, RUF=0.02, IR=2)
      
    # plotting above-ground conditions in minimum shade
      metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
    with(metout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                      , type = "l",main=paste("air and sky temperature",sep=""), ylim = c(-5, 40))})
    with(metout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                        , type = "l",lty=2,col='blue')})

    #vary parameters
    micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, Usrhyt=height, runshade = 0, spatial = spatial_path, minshade=25, maxshade=30, REFL=0.3, RUF=0.02, IR=2)
    
    # plotting above-ground conditions in minimum shade
    metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
    with(metout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                      , type = "l",col='darkorange',main=paste("air and sky temperature",sep=""), ylim = c(-5, 40))})
    with(metout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (Â°C)"
                        , type = "l",lty=2,col='cadetblue')})
    
      
    
  
  
  
  
  
  
