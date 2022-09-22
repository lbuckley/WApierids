library(ggplot2)
library(reshape2)
library(reshape)
library(viridisLite)
library(patchwork)
library(TrenchR)
library(tidyverse)

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

# stack models
d_stack <- dplyr::select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:weibull)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  dplyr::select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  dplyr::select(-fit) %>%
  unnest(preds)

# plot fits
ggplot(d_preds, aes(temp, rate)) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  facet_wrap(~model_name, labeller = labeller(model_name = label_facets_num), scales = 'free', ncol = 5) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none',
        strip.text = element_text(hjust = 0),
        strip.background = element_blank()) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Fits of every model available in rTPC') +
  geom_hline(aes(yintercept = 0), linetype = 2)

#extract model
mod= d_fits$beta[[1]]

## get predictions
#preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
#preds <- broom::augment(mod, newdata = preds)

#extract coefficients
tpc.beta= coef(mod)
tpc.beta$pop= pops[loc.k]

if(loc.k==1) tpc.betas= tpc.beta
if(loc.k>1) tpc.betas= rbind(tpc.betas, tpc.beta)

} #end location look

tpc.betas= as.data.frame(tpc.betas)

#plot tpcs
temps= 1:50

for(loc.k in 1:4){

ps= cbind(temps,
beta_2012(temps, as.numeric(tpc.betas[loc.k,1]), as.numeric(tpc.betas[loc.k,2]), as.numeric(tpc.betas[loc.k,3]), as.numeric(tpc.betas[loc.k,4]), as.numeric(tpc.betas[loc.k,5]) ))
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
#ASK FOR DATA
#analyze feeding rate distributions?

#plot temperature distributions
locations= c("Montrose","Sacramento", "LosBanos")

for(loc.k in c(1:2)){

#years for data
if (loc.k==1) years=c(1961:1971, 2001:2014, 2016:2021) 
#check 2015
if (loc.k==2) years=c(1961:1971, 2001:2016) 
if (loc.k==3) years=c(1962:1971, 2009:2018)

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

#--------------------------  
#plot density distributions
p1= ggplot(dat.day, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Feeding rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.5, 0.6),legend.background = element_rect(fill="transparent"))
#D0cm, TALOC, TAREF

#plot reference temperatures
p1.ref= ggplot(dat.day, aes(x=TAREF))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Feeding rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.5, 0.6),legend.background = element_rect(fill="transparent"))

if(loc.k==1) {temp.co= p1; temp.co.ref= p1.ref; dat.day.co=dat.day}
if(loc.k==2) {temp.ca= p1; temp.ca.ref= p1.ref; dat.day.ca=dat.day}

} # end location loop

#plot with TPCs
p1.all$TALOC= p1.all$temps

##FE paper fits
#co.colias= temp.co + geom_line(data=p1.all[p1.all$population=="CO",], aes(x=temps, y=p/5, lty=year),size=1.1)
##divide by 5 to line up, use shade temps? 
#ca.colias= temp.ca + geom_line(data=p1.all[p1.all$population=="CA",], aes(x=temps, y=p/5, lty=year),size=1.1)

#beta fits
co.colias= temp.co + geom_line(data=ps.all[ps.all$Population=="CO",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)
#divide by 5 to line up, use shade temps? 
ca.colias= temp.ca + geom_line(data=ps.all[ps.all$Population=="CA",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)

#add data
co.colias= co.colias + geom_point(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, y=mean/5, shape=Time),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5))
ca.colias= ca.colias + geom_point(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, y=mean/5, shape=Time),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5))

#reference temperatures
temp.co.ref + geom_line(data=ps.all[ps.all$Population=="CO",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)
temp.ca.ref + geom_line(data=ps.all[ps.all$Population=="CA",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)

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
topt.low= 25
topt.high= 45

if(loc.k==1) dat.day= dat.day.co
if(loc.k==2) dat.day= dat.day.ca

tpc.beta= unlist(tpc.betas[which(tpc.betas$pop==pop.times[loc.k]),1:5])

#estimate growth and development
#RUN
topts= seq(topt.low, topt.high, 1)
perf.mat= matrix(NA, nrow= nrow(dat.day), ncol= length(topts) )
for(topt.k in 1:length(topts)){
  perf.mat[,topt.k]= sapply(dat.day$TALOC, FUN=beta_2012, 
                            a=tpc.beta[1],b=topts[topt.k],c=tpc.beta[3],
                            d=tpc.beta[4],e=tpc.beta[5] )
  #set NaN to zero
  perf.mat[which(is.nan(perf.mat[,topt.k])),topt.k]= 0
}

#combine dates
perfs= cbind(dat.day[,c("year","period","seas")], perf.mat)

#aggregate
perfs1= aggregate(perfs[,4:24], list(perfs$year,perfs$period,perfs$seas), FUN=mean)
names(perfs1)=c("year","period","seas", seq(topt.low, topt.high, 1))

##not seasons, just study period
#perfs1= aggregate(perfs[,4:ncol(perfs)], list(perfs$year,perfs$period), FUN=mean)
#names(perfs1)=c("year","period",seq(topt.low, topt.high, 1))

#make column for slopes
perfs1$B= NA

for(row.k in 1: nrow(perfs1)){
  #put into format for regression
  perfs2= as.data.frame(cbind(topt.low:topt.high, t(perfs1[row.k,4:ncol(perfs1)]) ))
  colnames(perfs2)=c("topt","perf")
  
  #plot(perfs2$topt, perfs2$perf)
  #mod1= lm(perf~topt , data=perfs2)
  mod1= lm(perf~topt +I(topt^2) , data=perfs2)
  perfs1$B[row.k]= coefficients(mod1)[2]
  perfs1$curve[row.k]= coefficients(mod1)[3]
} 

#plot fitness curves through time
perfs.l= melt(perfs1[,1:24], id.vars = c("year","period","seas"))
names(perfs.l)[4:5]=c("temperature","performance")

#perfs.l= melt(perfs1[,1:ncol(perfs1)], id.vars = c("year","period"))
#names(perfs.l)[3:4]=c("temperature","performance")
perfs.l$temperature= as.numeric(as.character(perfs.l$temperature))

fig.fitnesscurves=ggplot(perfs.l, aes(x=temperature, y=performance, color=year, group=year))+geom_line()+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+
  xlab("thermal optima (C)")+ylab("growth rate (g/g/h)")+theme(legend.position = c(0.8, 0.3))

#find thermal optima through years 
perfs1$Topt= apply(perfs1[,4:ncol(perfs1)], MARGIN=1, FUN=which.max )

fig.topt= ggplot(perfs1, aes(x=year, y=Topt, color=seas))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)

if(loc.k==1) fig.fitnesscurves.co= fig.fitnesscurves
if(loc.k==2) fig.fitnesscurves.ca= fig.fitnesscurves

if(loc.k==1) fig.topt.co= fig.topt
if(loc.k==2) fig.topt.ca= fig.topt

} # end loop locations

#-------------------
#PLOT

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_Colias.pdf",height = 14, width = 14)
(co.colias+ ca.colias)/
  (fig.fitnesscurves.co + fig.fitnesscurves.ca)/
  (fig.topt.co + fig.topt.ca)+
  plot_annotation(tag_levels = 'A')
dev.off()

