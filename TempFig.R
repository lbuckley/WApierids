library(rnoaa)
library(trenchR)
library(lubridate)
library(reshape)
library(ggplot2)
library(dplyr)

#------------------------------------
#FIGURE 2: Temp

ghcnd(stationid = "AGE00147704")

#find stations: https://ncics.org/portfolio/monitor/ghcn-d-station-data/ 

# SEATTLE SAND POINT WEATHER FORECAST OFFICE, WA US (USW00094290)
# HATTON 9 SE, WA US (USC00453546)
# ALT: PRIEST RAPIDS DAM, WA US (USC00456747)

#western Colorado (Montrose and Mesa counties; elevation, 1,370-1,700 m
#grand junction USW00023066

sea= ghcnd_search("USW00094290", var = c("TMIN","TMAX"))
corfu= ghcnd_search("USC00453546", var = c("TMIN","TMAX"))
co= ghcnd_search("USW00023066", var = c("TMIN","TMAX"))
#ch= ghcnd_search("USW00013722", var = c("TMIN","TMAX"))  #"USC00311677" chapel hill, but issues with nrows

#assemble data
for(site in 1:3){
  
  #select time periods
  if(site==1) dat=sea
  if(site==2) dat=corfu
  if(site==3) dat=co
  
  t.dat=dat$tmin
  t.dat$tmin= t.dat$tmin/10 #divide by ten
  t.dat$tmax= dat$tmax$tmax/10
  t.dat$month= round(month(as.POSIXlt(t.dat$date)))
  t.dat$year= year(as.POSIXlt(t.dat$date))
  #restrict to growing season or summer
  #t.dat= t.dat[which(t.dat$month %in% c(4,5,7,8)),] 
  t.dat= t.dat[which(t.dat$month %in% c(7,8)),]
  
  #code season
  #month 3,4 vs 7,8
  t.dat$season<- NA
  t.dat$season[which(t.dat$month %in% c(4:5))] ="spring"
  t.dat$season[which(t.dat$month %in% c(7:8))] ="summer"
  
  #restrict years
  t.dat= t.dat[which(t.dat$year %in% c(1987:2021)),]
  t.dat1= t.dat[which(t.dat$year %in% c(1991:1994)),]
  #or 87-92, 1990:1994
  t.dat1$period="initial"
  t.dat2= t.dat[which(t.dat$year %in% c(2016,2018,2019,2020)),]
  t.dat2$period="recent"
  #combine
  t.dat= rbind(t.dat1,t.dat2)
  
  if(site==1){
    dat.sea=t.dat
  }
  if(site==2) {
    dat.cor=t.dat
  }
  if(site==3) {
    dat.co=t.dat
  }
}

dat.sea$site="Seattle, WA"
dat.cor$site="Corfu, WA"
dat.co$site="Western CO"
#dat.ch$site="Chapel Hill"
#combine
t.dat= rbind(dat.sea,dat.cor,dat.co)

#counts by year 
counts= t.dat %>% count(year, site)
#corfu: 1987:1994 complete, complete: 2016, 2018:2020
#sea: 1987:1989,1991:1994; 2014:2020

#scale to plant height
#don't scale since lacking surface temp?
air_temp_at_height_z<-function(z_0, z_r, z, T_r, T_s){
  T_r= as.numeric(T_r)
  T_s= as.numeric(T_s)
  T_z<-(T_r-T_s)*log((z+z_0)/z_0+1)/log((z_r+z_0)/z_0+1)+T_s ##this is exactly eqn (19) of the notes
  return(T_z)
}

#apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z_mat, z_0=z_0_1, z_r=2, z=z_0_1)}

#sine interpolation
#https://github.com/trenchproject/TrenchR/blob/master/R/DTRFunctions.R

#just day
dtr=function(T_max, T_min, t=7:18){
  gamma= 0.44 - 0.46* sin(0.9 + pi/12 * t)+ 0.11 * sin(0.9 + 2 * pi/12 * t);   # (2.2) diurnal temperature function
  T = T_max*gamma + T_min - T_min*gamma
  return(T)
}

temps= sapply(t.dat$tmax, FUN="dtr", T_min=t.dat$tmin)
temps= as.data.frame(t(temps))
temps$period= t.dat$period
temps$site= t.dat$site
temps$season= t.dat$season

#to long format
temps1<- melt(temps, id.vars=c("period","site","season"))

#change period names
temps1$period[temps1$period=="initial"]<-"1990s"
temps1$period[temps1$period=="recent"]<-"2010s"

#density plot
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/figures/")
pdf("TempPlot_78.pdf",height = 6, width = 6)
ggplot(temps1, aes(x=value, lty=period,color=site))+
  #facet_wrap(~season)+
  geom_density(lwd=1.2)+#scale_color_manual(values=c("orange","blue","darkgreen"))+
  xlim(0,45)+
  xlab("Temperature (°C)")+
  theme_classic(base_size = 18)+theme(legend.position = c(0.2, 0.7))+
  scale_color_viridis_d(begin=0, end=0.9)
dev.off()

#min, max distributions
ggplot(t.dat, aes(x=tmin, lty=period,color=site))+
  facet_wrap(~season)+
  geom_density()+scale_color_manual(values=c("blue","orange","darkgreen"))+
  xlim(-20,30)+
  xlab("Temperature (°C)")+
  theme_bw(base_size = 18)

ggplot(t.dat, aes(x=tmax, lty=period,color=site))+
  facet_wrap(~season)+
  geom_density()+scale_color_manual(values=c("blue","orange","darkgreen"))+
  xlim(0,45)+
  xlab("Temperature (°C)")+
  theme_bw(base_size = 18)

#recent soil temperatures
#https://weather.wsu.edu/
 
 
