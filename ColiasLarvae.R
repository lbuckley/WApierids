library(ggplot2)
library(reshape2)
library(reshape)
library(viridisLite)
library(patchwork)
library(TrenchR)
library(tidyverse)
library(rTPC)
library(nls.multstart)
library(dplyr)

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

#add two weeks to account for lag between pupae and adults
oldData$doy= oldData$doy +14
contempSummary$doy= contempSummary$doy +14

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_ColiasPhoto.pdf",height = 6, width = 6)

#plot export dim (paper): 750x500 pixels
plot(oldData$doy,oldData$Reflectance,ylim=c(0.2,0.5),cex=1.5, cex.axis=1.2,cex.lab=1.3,ylab="Reflectance (at 650 nm)",xlab="Day of Year", bty="l")
arrows(oldData$doy, oldData$Reflectance-oldData$RefSE*1.96, oldData$doy, oldData$Reflectance+oldData$RefSE*1.96, length=0.05, angle=90, code=3,lwd=2)
lines(oldData$doy,oldData$Reflectance,lwd=2, lty=2) #range limits range of doys graphed

points(contempSummary$doy,contempSummary$means,col="blue",cex=1.5, pch=2)
arrows(contempSummary$doy, contempSummary$means-contempSummary$se*1.96, contempSummary$doy, contempSummary$means+contempSummary$se*1.96, length=0.05, angle=90, code=3,lwd=2, col="blue")
lines(contempSummary$doy,contempSummary$means,col="blue",lwd=2)
legend(100,.3,c(1971,2018),col=c("black","blue"),lwd=2,lty=c(2,1)) #presentation only

dev.off()

#dates
# c(0,50,100,150)
# Jan 1, Feb 19, Apr 10, May 30

#plot temperature distributions
locations= c("Sacramento")
loc.k=1

#===============================================================
#CO larvae

temps=0:50

#Colias butterfly larval feeding
#Higgins
tpc= function(T, Fmax,To, row, sigma) Fmax*exp(-exp(row*(T-To)-6)-sigma*(T-To)^2)

setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/Higgins/')
dat= read.csv("HigginsTPC.csv")
#keep comparison
dat= dat[c(1,3,5,6),]

for(i in 1:nrow(dat)){
  p= tpc(T=temps, Fmax=dat$Fmax[i],To=dat$Topt[i], row=dat$row[i], sigma=dat$sigma[i])
  p1=as.data.frame(cbind(temps,p))
  p1$taxa="butterfly"
  p1$species=dat$species[i]
  p1$population=dat$population[i]
  p1$year= dat$year[i]
  
  if(i==1) p1.all= p1
  if(i>1) p1.all=rbind(p1, p1.all)
}

#change population labels
p1.all$population[p1.all$population=="Sacramento Valley, CA"]="CA"
p1.all$population[p1.all$population=="Montrose Valley, CO"]="CO"

#plot
p1.all$year= as.factor(p1.all$year)

fig2.butterfly= ggplot(p1.all,aes(x=temps, y=p))+geom_line(aes(color=population, lty=year),size=1.1)+
  theme_bw(base_size=14)+xlab("")+ylab("")+
  theme(legend.position = c(0.25, 0.65),legend.background = element_rect(fill="transparent"))+
  scale_color_viridis_d()

#------------------------
#plot from Higgins data
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/Higgins/')
rdat= read.csv("totalFT.csv")
hdat= read.csv("FR1972.csv")

se <- function(x) sqrt(var(x)/length(x))

#historic
#all
ggplot(hdat,aes(x=Temp, y=FR, color=Pop))+geom_point() +geom_smooth()

#temp means and se
hdat.m= hdat %>%                                        
  group_by(Pop, Temp) %>%                        
  summarise_at(vars(FR),
               list(mean = mean, se=se))

ggplot(hdat.m,aes(x=Temp, y=mean, color=Pop))+geom_point() +geom_smooth(se=F)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))

#----------
#current
#all
ggplot(rdat,aes(x=Temp, y=lnTMG, color=pop))+geom_point() +geom_smooth()

#temp means and se
rdat.m= rdat %>%                                        
  group_by(pop, Temp) %>%                        
  summarise_at(vars(lnTMG),
               list(mean = mean, se=se))

ggplot(rdat.m,aes(x=Temp, y=mean, color=pop))+geom_point() +geom_smooth(se=F)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se))

#----------
#Plot beta fits

colnames(rdat.m)[1]= "Pop"
hdat.m$Time= "historic"
rdat.m$Time= "recent"
#change Olathe name
rdat.m$Pop[rdat.m$Pop=="OL"]="CO"
#drop Gun and NC populations
rdat.m= rdat.m[rdat.m$Pop %in% c("CA","CO"),]

c.dat= rbind(hdat.m,rdat.m)
c.dat$PopTime= paste(c.dat$Pop, c.dat$Time, sep="_") 
pops= unique(c.dat$PopTime)

for(loc.k in 1:length(pops)){

d= c.dat[which(c.dat$PopTime==pops[loc.k]),c("Temp","mean")]
colnames(d)=c("temp","rate")

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(beta = map(data, ~nls_multstart(rate~beta_2012(temp = temp, a, b, c, d, e),
                                         data = .x,
                                         iter = c(6,6,6,6,6),
                                         start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') - 10,
                                         start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'beta_2012') + 10,
                                         lower = get_lower_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         upper = get_upper_lims(.x$temp, .x$rate, model_name = 'beta_2012'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         weibull = map(data, ~nls_multstart(rate~weibull_1995(temp = temp, a,topt,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'weibull_1995') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'weibull_1995'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)))
} #end loc.k loop

#-------------
#Fit Asbury and Angilletta beta function

for(loc.k in 1:length(pops)){
  
  d= c.dat[which(c.dat$PopTime==pops[loc.k]),c("Temp","mean")]
  colnames(d)=c("temp","rate")

#scale to max height
d$rate_scale= d$rate*2.5/max(d$rate)

fits <- nls(rate_scale~TPC_beta(T_b=temp, shift, breadth, tolerance, aran=0, skew),
            data = d,
            start= c(shift = 0, breadth = 0.1, tolerance=43, skew=0.6), 
            algorithm = "port",
            lower = c(shift = -9, breadth = 0.04, tolerance=20, skew=0.5),
            upper = c(shift = 9, breadth = 0.15, tolerance=60, skew=0.8)
        )

# param.grid= expand.grid(shift=-9:5, breadth=seq(0.04,0.14,0.01), skew=seq(0.1,0.9,0.1) )
# plot(1:50, TPC_beta(T_b=1:50, shift=-9, breadth=0.07, aran=0, tolerance=43, skew=0.7), type="l")
# for(n.row in 1:nrow(param.grid)){
#   points(1:50, TPC_beta(T_b=1:50, shift=param.grid[n.row,1], breadth=param.grid[n.row,2], aran=0, tolerance=43, skew=0.7), type="l")
# }

coefs= summary(fits)$coefficients

#extract coefficients
tpc.beta= coef(fits)
tpc.beta$pop= pops[loc.k]
tpc.beta$maxp= max(d$rate)

if(loc.k==1) tpc.betas= tpc.beta
if(loc.k>1) tpc.betas= rbind(tpc.betas, tpc.beta)

} #end location loop

tpc.betas= as.data.frame(tpc.betas)
tpc.betas$shift= as.numeric(tpc.betas$shift)
tpc.betas$breadth= as.numeric(tpc.betas$breadth)
tpc.betas$tolerance= as.numeric(tpc.betas$tolerance)
tpc.betas$skew= as.numeric(tpc.betas$skew)
tpc.betas$maxp= as.numeric(tpc.betas$maxp)

#plot with points
par(mfrow=c(2,2))

for(loc.k in 1:4){
  d= c.dat[which(c.dat$PopTime==pops[loc.k]),c("Temp","mean")]
  colnames(d)=c("temp","rate")
  
  #scale to max height
  d$rate_scale= d$rate*2.5/max(d$rate)
  
  plot(d$temp,d$rate_scale, ylim=c(0,3), xlim=c(0,50))
  points(1:50, TPC_beta(T_b=1:50, shift=tpc.betas[loc.k,1], breadth=tpc.betas[loc.k,2], aran=0, tolerance=tpc.betas[loc.k,3], skew=tpc.betas[loc.k,4]), type="l")
  
}

#-------------

#plot tpcs
temps= 1:50

for(loc.k in 1:4){

ps= cbind(temps,
          TPC_beta(T_b=1:50, shift=tpc.betas[loc.k,1], breadth=tpc.betas[loc.k,2], aran=0, tolerance=tpc.betas[loc.k,3], skew=tpc.betas[loc.k,4])/(2.5/tpc.betas[loc.k,6]) )
ps= as.data.frame(ps)
ps$Pop= pops[loc.k]

if(loc.k==1) ps.all=ps
if(loc.k>1) ps.all= rbind(ps.all, ps)
}

colnames(ps.all)= c("Temp","Perf","Pop")
ps.all$Perf[is.nan(ps.all$Perf)]=0

ps.all$Population= "CO"
ps.all$Population[grep("CA", ps.all$Pop)]="CA"
ps.all$Time= "historic"
ps.all$Time[grep("recent", ps.all$Pop)]="recent"

#plot
ggplot(ps.all,aes(x=Temp, y=Perf, color=Pop))+geom_line()

#==============================================================
#Colias feeding

#Fig 1: modern
#Fig 3: Historic, Sacramento and Montrose; 1961-1971, 2001-2011 
#analyze feeding rate distributions?

#plot temperature distributions
locations= c("Montrose","Sacramento", "LosBanos")

for(loc.k in c(1:2)){

#years for data
if (loc.k==1) years=c(1961:2014, 2016:2021) 
#check 2015
if (loc.k==2) years=c(1961:2021) 
if (loc.k==3) years=c(1962:1971, 2009:2018)

#SUN
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/')
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

#drop middle period
dat.day.plot= dat.day[is.na(dat.day$period),]

#--------------------------  
#plot density distributions
p1= ggplot(dat.day.plot, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Feeding rate (g/g/h)")+
  xlab("Temperature (°C)" )+
  xlim(-5,50)+
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.082), expand=c(0,0))
#D0cm, TALOC, TAREF

#plot reference temperatures
p1.ref= ggplot(dat.day.plot, aes(x=TAREF))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Feeding rate (g/g/h)")+
  xlab("Temperature at reference height (°C)" )+
  xlim(-5,50)+ 
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.6, 0.7),legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(limits=c(0,0.082), expand=c(0,0))

#plot together
p1.pl.ref=p1 + geom_density(alpha=0.4, linetype="dashed", aes(x=TAREF, color=period))+
  theme(legend.position = c(0.6, 0.6))

#remove label
p1= p1+theme(strip.text.x = element_blank())

if(loc.k==1) {temp.co= p1; temp.co.ref= p1.ref; dat.day.co=dat.day; temps.co= p1.pl.ref}
if(loc.k==2) {temp.ca= p1; temp.ca.ref= p1.ref; dat.day.ca=dat.day; temps.ca= p1.pl.ref}

} # end location loop

#plot with TPCs
p1.all$TALOC= p1.all$temps

##FE paper fits
#co.colias= temp.co + geom_line(data=p1.all[p1.all$population=="CO",], aes(x=temps, y=p/5, lty=year),size=1.1)
##divide by 5 to line up, use shade temps? 
#ca.colias= temp.ca + geom_line(data=p1.all[p1.all$population=="CA",], aes(x=temps, y=p/5, lty=year),size=1.1)

#change names
ps.all$period= ps.all$Time
ps.all$period[ps.all$period=="historic"]= "1961-1971"
ps.all$period[ps.all$period=="recent"]= "2001-2011"
c.dat$period= c.dat$Time
c.dat$period[c.dat$period=="historic"]= "1961-1971"
ps.all$period[ps.all$period=="recent"]= "2001-2011"

#beta fits
co.colias= temp.co + geom_line(data=ps.all[ps.all$Population=="CO",], aes(x=Temp, y=Perf/5, color=period),size=1.1)
#divide by 5 to line up, use shade temps? 
ca.colias= temp.ca + geom_line(data=ps.all[ps.all$Population=="CA",], aes(x=Temp, y=Perf/5, color=period),size=1.1)

#add data
co.colias= co.colias + geom_point(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, y=mean/5, color=period),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5, color=period))
ca.colias= ca.colias + geom_point(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, y=mean/5, color=period),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5, color=period))

#reference temperatures
co.colias.ref= temp.co.ref + geom_line(data=ps.all[ps.all$Population=="CO",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)+
  geom_point(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, y=mean/5, shape=Time),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5))
ca.colias.ref= temp.ca.ref + geom_line(data=ps.all[ps.all$Population=="CA",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)+
  geom_point(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, y=mean/5, shape=Time),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5))

#--------------------
#Global hourly 
#https://www.ncei.noaa.gov/maps/hourly/
# grand junction: 72476023066
# travis field: 74516023202

#compare to GHCN station data
#Winters, CA: GHCND:USC00049742
#Montrose 2, CO: GHCND:USC00055722

#library(rnoaa)
#wint<- ghcnd(stationid = "USC00049742")
#wint= filter(wint, element %in% c("TMIN","TMAX"), year %in% c(1961:1971, 2001:2021) )

#-------------------
#CO selection analysis

#a= height; b= mode; c= breadth; d= inverse breadth?; e= breadth

pop.times= c("CO_historic","CA_historic")

for(loc.k in 1:length(pop.times)){

#generate parameter combinations
tpc.beta= as.numeric(tpc.betas[which(tpc.betas$pop==pop.times[loc.k]),c(1:4,6)] ) 

param.grid= expand.grid(shift= seq(tpc.beta[1]-10,tpc.beta[1]+10,1),
                        breadth=c(tpc.beta[2]-0.03,tpc.beta[2],tpc.beta[2]+0.03) )

if(loc.k==1) dat.day= dat.day.co
if(loc.k==2) dat.day= dat.day.ca

#estimate growth and development
#RUN
perf.mat= matrix(NA, nrow= nrow(dat.day), ncol= nrow(param.grid) )
for(topt.k in 1:nrow(param.grid) ){
  perf.mat[,topt.k]= sapply(dat.day$TALOC, FUN=TPC_beta, 
                            shift=param.grid[topt.k,1], breadth=param.grid[topt.k,2], aran=0, 
                            tolerance=tpc.beta[3], skew=tpc.beta[4])
  
  #set NaN to zero
  perf.mat[which(is.nan(perf.mat[,topt.k])),topt.k]= 0
}

#combine dates
perfs= cbind(dat.day[,c("year","period","seas")], perf.mat)

#aggregate
perfs1= aggregate(perfs[,4:ncol(perfs)], list(perfs$year,perfs$period,perfs$seas), FUN=mean)
names(perfs1)=c("year","period","seas", 1:nrow(param.grid) )

##not seasons, just study period
#perfs1= aggregate(perfs[,4:ncol(perfs)], list(perfs$year,perfs$period), FUN=mean)
#names(perfs1)=c("year","period",seq(topt.low, topt.high, 1))

# #make column for slopes
# perfs1$B= NA
# 
# for(row.k in 1: nrow(perfs1)){
#   #put into format for regression
#   perfs2= as.data.frame(cbind(topt.low:topt.high, t(perfs1[row.k,4:ncol(perfs1)]) ))
#   colnames(perfs2)=c("topt","perf")
#   
#   #plot(perfs2$topt, perfs2$perf)
#   #mod1= lm(perf~topt , data=perfs2)
#   mod1= lm(perf~topt +I(topt^2) , data=perfs2)
#   perfs1$B[row.k]= coefficients(mod1)[2]
#   perfs1$curve[row.k]= coefficients(mod1)[3]
# } 

#estimate topt for curve
param.grid$Topt=NA
for(topt.k in 1:nrow(param.grid)){
  ps= TPC_beta(1:50,shift=param.grid[topt.k,1], breadth=param.grid[topt.k,2], aran=0, 
               tolerance=tpc.beta[3], skew=tpc.beta[4])
  param.grid$Topt[topt.k]= c(1:50)[which.max(ps)]
}

#plot fitness curves through time
perfs.l= melt(perfs1[,1:(ncol(perfs)-1)], id.vars = c("year","period","seas"))
names(perfs.l)[4:5]=c("temperature","performance")

perfs.l$shift= param.grid$shift[perfs.l$temperature]
perfs.l$breadth= param.grid$breadth[perfs.l$temperature]
perfs.l$topt= param.grid$Topt[perfs.l$temperature]

#perfs.l= melt(perfs1[,1:ncol(perfs1)], id.vars = c("year","period"))
#names(perfs.l)[3:4]=c("temperature","performance")
perfs.l$temperature= as.numeric(as.character(perfs.l$temperature))

#group by breadth
perfs.l$yrbr= paste(perfs.l$year, perfs.l$breadth,sep="_")
perfs.l$breadth= factor(perfs.l$breadth)
#scale performance
perfs.l$performance= perfs.l$performance/(2.5/tpc.beta[5])

#combine shift, breadth
fig.fitnesscurves=ggplot(perfs.l, aes(x=topt, y=performance, color=year, group=yrbr) )+geom_line(aes(lty=breadth))+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+ xlim(-5,50)+
  xlab("Thermal optima (C)")+ylab("Mean feeding rate (g/g/h)") #+theme(legend.position = c(0.8, 0.3))

#fitness at fixed breadth
fig.fit.fb=ggplot(perfs.l[perfs.l$breadth==0.15,], aes(x=topt, y=performance, color=year, group=year) )+geom_line()+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+ xlim(-5,50)+
  xlab("Thermal optima (C)")+ylab("Mean feeding rate (g/g/h)")+theme(legend.position = c(0.6, 0.8))+
  theme(strip.text.x = element_blank())

#fitness at fixed shift
fig.fit.fs=ggplot(perfs.l[round(perfs.l$shift,3)== 3.652,], aes(x=breadth, y=performance, color=year, group=year) )+geom_line()+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+ xlim(-5,50)+
  xlab("Thermal optima (C)")+ylab("Mean feeding rate (g/g/h)")+theme(legend.position = c(0.8, 0.3))

#find thermal optima through years  
breadths= sort(unique(param.grid$breadth))

inds= which(param.grid$breadth==breadths[1])+3
perfs1$shift_opt_b1= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

inds= which(param.grid$breadth==breadths[2])+3
perfs1$shift_opt_b2= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

inds= which(param.grid$breadth==breadths[3])+3
perfs1$shift_opt_b3= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

#melt 
perfs.b= perfs1[,c("year","period","seas","shift_opt_b1","shift_opt_b2","shift_opt_b3")]
perfs.b.l= melt(perfs.b, id.vars = c("year","period","seas"))
names(perfs.b.l)[4:5]=c("breadth","opt_shift")
perfs.b.l$breadths= NA
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b1"]<- breadths[1]
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b2"]<- breadths[2]
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b3"]<- breadths[3]

perfs.b.l$seas_br= paste(perfs.b.l$seas, perfs.b.l$breadths,sep="_")
perfs.b.l$breadth= factor(perfs.b.l$breadths)

#beta= 0.15
fig.shift_opt= ggplot(perfs.b.l[which(perfs.b.l$breadth==0.15),], aes(x=year, y=opt_shift, color=seas, group= seas_br, lty=breadth))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  xlab("Year")+ylab("Optimal thermal optima (C)")+
  guides(lty="none")+theme(legend.position = c(0.9, 0.7))+ 
  scale_color_manual(values= c("orange","purple"))+
  labs(colour = "season") 

#all betas
fig.shift_opt.all= ggplot(perfs.b.l, aes(x=year, y=opt_shift, color=seas, group= seas_br, lty=breadth))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  xlab("Year")+ylab("Optimal thermal optima (C)")+
  scale_color_manual(values= c("orange","purple"))+
  labs(colour = "season") 

if(loc.k==1) {fig.fitnesscurves.co= fig.fit.fb; fig.fit.all.co= fig.fitnesscurves}
if(loc.k==2) {fig.fitnesscurves.ca= fig.fit.fb; fig.fit.all.ca= fig.fitnesscurves}

if(loc.k==1) {fig.shift_opt.co= fig.shift_opt; fig.shift_opt.all.co= fig.shift_opt.all}
if(loc.k==2) {fig.shift_opt.ca= fig.shift_opt; fig.shift_opt.all.ca= fig.shift_opt.all}

} # end loop locations

#-------------------
#PLOT

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_Colias_CO.pdf",height = 18, width = 8)
temps.co/
  co.colias/
  fig.fitnesscurves.co/
  fig.shift_opt.co +
  plot_annotation(tag_levels = 'A')
dev.off()

pdf("Fig_Colias_CA.pdf",height = 18, width = 8)
temps.ca/
  ca.colias/
  fig.fitnesscurves.ca/
  fig.shift_opt.ca +
  plot_annotation(tag_levels = 'A')
dev.off()

#supplementary plots
pdf("Fig_Colias_CO_supp.pdf",height = 10, width = 8)
fig.fit.all.co/
  fig.shift_opt.all.co +
  plot_annotation(tag_levels = 'A')
dev.off()

pdf("Fig_Colias_CA_supp.pdf",height = 10, width = 8)
fig.fit.all.ca/
  fig.shift_opt.all.ca +
  plot_annotation(tag_levels = 'A')
dev.off()




