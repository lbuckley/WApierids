library(ggplot2)
library(reshape2)
library(reshape)
library(viridisLite)
library(patchwork)
library(TrenchR)

#Montrose Valley, CO: 1961-1971; 2001-2011;N 38.62, W 108.02, 1633m
#Sacramento Valley, CA: 1961-1971; 2001-2011; N 38.44, W121.86, 19m

#Empirical data
#Nielsen and Kingsolver
#https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13515
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/NielsenColias/')

oldData= read.csv("Nielsen_oldData.csv")
#cut data to matching photoperiods
oldData= oldData[which(oldData$Photoperiod %in% 10:16),]

contempSummary= read.csv("Nielsen_contempSummary.csv")

#convert photoperiod to days of year
daylengths= daylength(lat=37.983, doy=1:173)

a= oldData$Photoperiod
b= daylengths
oldData$doy= sapply(a, function(a, b) {which.min(abs(a-b))}, b)

a= contempSummary$Larval_Photoperiod
b= daylengths
contempSummary$doy= sapply(a, function(a, b) {which.min(abs(a-b))}, b)

#plot export dim (paper): 750x500 pixels
plot(oldData$doy,oldData$Reflectance,ylim=c(0.2,0.5),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Reflectance (at 650 nm)",xlab="Day of Year", bty="l")
arrows(oldData$doy, oldData$Reflectance-oldData$RefSE*1.96, oldData$doy, oldData$Reflectance+oldData$RefSE*1.96, length=0.05, angle=90, code=3,lwd=2)
lines(oldData$doy,oldData$Reflectance,lwd=2, lty=2) #range limits range of doys graphed

points(contempSummary$doy,contempSummary$means,col="blue",cex=1.5, pch=2)
arrows(contempSummary$doy, contempSummary$means-contempSummary$se*1.96, contempSummary$doy, contempSummary$means+contempSummary$se*1.96, length=0.05, angle=90, code=3,lwd=2, col="blue")
lines(contempSummary$doy,contempSummary$means,col="blue",lwd=2)
legend(100,.3,c(1971,2018),col=c("black","blue"),lwd=2,lty=c(2,1)) #presentation only

#dates
# c(0,50,100,150)
# Jan 1, Feb 19, Apr 10, May 30

#plot temperature distributions
locations= c("Sacramento")
loc.k=1

#==============================================================
#Montrose Colias

#Fig 1: modern
#Fig 3: Historic, Sacramento and Montrose; 1961-1971, 2001-2011
#ASK FOR DATA
#analyze feeding rate distributions?

#plot temperature distributions
loc.k=2
locations= c("Montrose","Sacramento")

#years for data
if (loc.k==1) years=c(1961:1971, 2001:2014, 2016:2021) 
#check 2015
if (loc.k==2) years=c(1961:1971, 2001:2016) 

#SUN
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/')
#combine data
for(yr.k in 1:length(years)){
dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
dat$year= years[yr.k]

#vary periods
  if(years[yr.k]<2001)dat$period="initial"
  if(years[yr.k]>=2001 & years[yr.k]<2011)dat$period="middle"
  if(years[yr.k]>=2011)dat$period="recent"

if(yr.k==1) dat.all=dat
if(yr.k>1) dat.all=rbind(dat, dat.all)
}

#subset to sunlight
dat.day= subset(dat.all, dat.all$SOLR>0)

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

dat.day.sun= dat.day

#--------------------------  
#plot density distributions
p1= ggplot(dat.day.sun, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Growth rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.5, 0.8))
#D0cm, TALOC, TAREF

