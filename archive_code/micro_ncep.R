# Adapt Colias niche model

#global historical climatology network?

#Center for urban horticulture
#lat 47.657628, Longitude: -122.290255

#Corfu site in central Washington in the Columbia National Wildlife Refuge
#12 miles west of Othello
#46.850663264 -119.535331192
  
#------------------

library(NicheMapR)
library(RNCEP)
library(elevatr)
library(microclima)

micro <- micro_ncep(loc = lonlat, dstart = dstart, dfinish = dfinish, DEP = DEP,
                    runmoist = 0, runshade = 0, Usrhyt = 0.01)

metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade

# append dates
tzone<-paste("Etc/GMT+",0,sep="")
dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")

metout <- cbind(dates,metout)
soil <- cbind(dates,soil)
soilmoist <- cbind(dates, soilmoist)

# plotting above-ground conditions in minimum shade
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

#extract air temperature and soil
dat <- cbind(metout[,1:3],metout$TAREF, metout$ZEN, metout$SOLR, soil$D0cm)


