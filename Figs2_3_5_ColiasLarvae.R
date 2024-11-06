#R version 4.1.0 (2021-05-18)

library(ggplot2) #ggplot2_3.4.4
library(reshape2) #reshape2_1.4.4
library(reshape) #reshape_0.8.8
library(viridisLite) #viridisLite_0.4.0
library(patchwork) #patchwork_1.2.0
library(TrenchR) #TrenchR_1.1.1
library(tidyverse) #tidyverse_1.3.2
library(rTPC) #rTPC_1.0.2
library(nls.multstart) #nls.multstart_1.2.0
library(dplyr) #dplyr_1.1.2

#Montrose Valley, CO: 1961-1971; 2001-2011;N 38.62, W 108.02, 1633m
#Sacramento Valley, CA: 1961-1971; 2001-2011; N 38.44, W121.86, 19m

#Empirical data
#Nielsen and Kingsolver
#https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13515

oldData= read.csv("./data/NielsenColias/Nielsen_oldData.csv")
#cut data to matching photoperiods
oldData= oldData[which(oldData$Photoperiod %in% 10:16),]

contempSummary= read.csv("./data/NielsenColias/Nielsen_contempSummary.csv")

#convert photoperiod to days of year
daylengths= daylength(lat=37.983, doy=1:173)

a= oldData$Photoperiod
b= daylengths
oldData$doy= sapply(a, function(a, b) {which.min(abs(a-b))}, b)

a= contempSummary$Larval_Photoperiod
b= daylengths
contempSummary$doy= sapply(a, function(a, b) {which.min(abs(a-b))}, b)

#Use development data to account for lag between pupae and adults
#load development data
dt=  read.csv("./data/NielsenColias/DevelopmentTime_Table2_AllenSmith.csv")

#use temperature estimates to estimate pupation time
plot(dt$Temperature, 1/dt$DT_pupa)
dt.mod= lm(1/dt$DT_pupa ~dt$Temperature)
b= dt.mod$coefficients[1]
a= dt.mod$coefficients[2]
To= -b/a
DD= 1/a

#extract temperature data
locations= c("Montrose","Sacramento", "LosBanos")
yearsn= c(1970, 2018)
doys= oldData$doy

dat1970= read.csv("./data/era5_micro_sun/Sacramento1970_Jan.csv")
dat2018= read.csv("./data/era5_micro_sun/Sacramento2018_Jan.csv")

dat1970$dd= dat1970$TAREF -To
dat2018$dd= dat2018$TAREF -To
dat1970$dd[dat1970$dd<0]=0
dat2018$dd[dat2018$dd<0]=0

dat1970$cdd= cumsum(dat1970$dd)/24
dat2018$cdd= cumsum(dat2018$dd)/24

for(yr.k in 1:2){
  doy.adult= rep(NA, length(doys))
  
  if(yr.k==1) dat=dat1970
  if(yr.k==2) dat=dat2018
  
for(day.k in 1:length(doys)){
  
  dat$cdd.doy= dat$cdd - dat[match(doys[day.k],dat$DOY),"cdd"]
  doy.adult[day.k] = dat[which.max(dat$cdd.doy>DD), "DOY"]
}
  if(yr.k==1) doy.adult.1970= doy.adult #Fix to 1971
  if(yr.k==2) doy.adult.2018= doy.adult
}

#add days
match1= match(oldData$doy, doys)
oldData$AdultDOY= doy.adult.1970[match1]

match1= match(contempSummary$doy, doys)
contempSummary$AdultDOY= doy.adult.2018[match1]

#check development time
#cap pupal time at 40 days
old_Dt=oldData$AdultDOY - oldData$doy
cont_dt=contempSummary$AdultDOY - contempSummary$doy
oldData$AdultDOY[which(old_Dt>40)]= oldData$doy[which(old_Dt>40)]+40
contempSummary$AdultDOY[which(cont_dt>40)]= oldData$doy[which(cont_dt>40)]+40

#---

#plot ggplot
od= oldData[,c("doy","Reflectance", "RefSE","AdultDOY")]
cd= contempSummary[,c("doy","means", "se","AdultDOY")]
names(cd)= names(od)
od$year= 1971
cd$year= 2018
od= rbind(od, cd)
od$year= factor(od$year)
od$AdultDOY1971= c(oldData$AdultDOY, oldData$AdultDOY[1:5])

#Fig 5
fig.co.photo= ggplot(od,aes(x=AdultDOY1971, y=Reflectance, color=year, group=year))+ geom_segment(aes(x = 60, y =0.48, xend =120, yend =0.48))+
  geom_segment(aes(x = 182, y =0.48, xend =185, yend =0.48), arrow = arrow(length = unit(0.1, "inches")), show.legend = FALSE)+
  annotate("text", x = 90, y = 0.49, label = "Mar Apr")+
  annotate("text", x = 175, y = 0.49, label = "Jul Aug")

fig.co.photo= fig.co.photo+geom_line(linewidth=1)+geom_point()+
  geom_errorbar(aes(x=AdultDOY1971, ymin=Reflectance-RefSE, ymax=Reflectance+RefSE) )+
  theme_classic(base_size = 20)+
  theme(legend.position = c(0.7, 0.4),legend.background = element_rect(fill="transparent"))+
  scale_color_manual(values=c("#7AD151FF","#440154FF"))+
  xlab("Day of Year") +ylab("Reflectance (at 650 nm)")

#dates
# c(0,50,100,150)
# Jan 1, Feb 19, Apr 10, May 30

#===============================================================
#CO larvae

temps=0:50

#Colias butterfly larval feeding
#Data from Higgins JK, MacLean HJ, Buckley LB, Kingsolver JG. Geographic differences and microevolutionary changes in thermal sensitivity of butterfly larvae in response to climate. Functional Ecology. 2014 Aug 1:982-9. https://www.jstor.org/stable/24033587

tpc= function(T, Fmax,To, row, sigma) Fmax*exp(-exp(row*(T-To)-6)-sigma*(T-To)^2)

dat= read.csv("./data/Higgins/HigginsTPC.csv")
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
rdat= read.csv("./data/Higgins/totalFT.csv")
hdat= read.csv("./data/Higgins/FR1972.csv")

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

se<- coefs[,2]

if(loc.k==1) {tpc.betas= tpc.beta; tpc.se= se}
if(loc.k>1) {tpc.betas= rbind(tpc.betas, tpc.beta); tpc.se= rbind(tpc.se, se)}

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
ps.all$Temp<- as.numeric(ps.all$Temp)
tpc.plot<- ggplot(ps.all,aes(x=Temp, y=Perf, color=Pop))+geom_line()

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
if (loc.k==2) years=c(1961:2016,2018:2021) 
if (loc.k==3) years=c(1962:1971, 2009:2018)

#combine data
for(yr.k in 1:length(years)){
dat= read.csv(paste("./data/era5_micro_sun/",locations[loc.k],years[yr.k],".csv",sep="") )

dat$year= years[yr.k]

dat$period<- NA

  if(years[yr.k]<1972) dat$period<- "1961-1971"
  if(years[yr.k]>=2001 & years[yr.k]<=2011) dat$period <- "2001-2011"

if(yr.k==1) dat.all=dat
if(yr.k>1) dat.all=rbind(dat, dat.all)
}

#keep all hours, no longer subset to sunlight
dat.day=dat.all 

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
  #ylab("Feeding rate (g/g/h)")+
  ylab("Density" )+
  xlab("Temperature (°C)" )+
  xlim(-5,50)+
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.075), expand=c(0,0))+
  scale_color_manual(values=c("#7AD151FF","#440154FF","black")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF","black")  ) 
#D0cm, TALOC, TAREF

#plot reference temperatures
p1.ref= ggplot(dat.day.plot, aes(x=TAREF))+
  geom_density(alpha=0.3, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  #ylab("Feeding rate (g/g/h)")+
  ylab("Density" )+
  xlab("Temperature at reference height (°C)" )+
  xlim(-5,50)+ 
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.4, 0.8),legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(limits=c(0,0.075), expand=c(0,0))+
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 

#plot together
p1.pl.ref=p1 + geom_density(alpha=0.3, linetype="dashed", aes(x=TAREF, color=period))+
  theme(legend.position = c(0.4, 0.8))

#remove label
p1= p1+theme(strip.text.x = element_blank())

# plot for Nielsen study
if(loc.k==2){
  dat.day$period[dat.day$year %in% c(2011:2021)] <- "2011-2021"
  
  #drop middle period
  dat.day.plot= dat.day[which(dat.day$period %in% c("1961-1971", "2011-2021")),]
  #make character factor
  dat.day.plot$period= as.factor(as.character(dat.day.plot$period))
  
  #--------------------------  
  #plot density distributions
  p1n= ggplot(dat.day.plot, aes(x=TALOC))+
    geom_density(alpha=0.3, aes(fill=period, color=period))+
    facet_wrap(~seas)+
    ylab("Density")+ #"Feeding rate (g/g/h)"
    xlab("Temperature (°C)" )+
    xlim(-5,50)+
    theme_classic(base_size = 20) +
    theme(legend.position = "none")+
    scale_y_continuous(limits=c(0,0.075), expand=c(0,0))+
    scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
    scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 
  
  #D0cm, TALOC, TAREF
  
  #plot reference temperatures
  p1n.ref= ggplot(dat.day.plot, aes(x=TAREF))+
    geom_density(alpha=0.3, aes(fill=period, color=period))+
    facet_wrap(~seas)+
    ylab("Density")+
    xlab("Temperature at reference height (°C)" )+
    xlim(-5,50)+ 
    theme_classic(base_size = 20) +
    theme(legend.position = c(0.6, 0.7),legend.background = element_rect(fill="transparent"))+
    scale_y_continuous(limits=c(0,0.075), expand=c(0,0))+
    scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
    scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 
  
  #plot together
  p1n.pl.ref=p1n + geom_density(alpha=0.3, linetype="dashed", aes(x=TAREF, color=period))+
    theme(legend.position = c(0.6, 0.85),legend.background = element_rect(fill="transparent")) 
}

if(loc.k==1) {temp.co= p1; temp.co.ref= p1.ref; dat.day.co=dat.day; temps.co= p1.pl.ref}
if(loc.k==2) {temp.ca= p1; temp.ca.ref= p1.ref; dat.day.ca=dat.day; temps.ca= p1.pl.ref; temp.nielsen=p1n.pl.ref}

} # end location loop

#plot with TPCs
p1.all$TALOC= p1.all$temps

#change names
ps.all$period= ps.all$Time
ps.all$period[ps.all$period=="historic"]= "1961-1971"
ps.all$period[ps.all$period=="recent"]= "2001-2011"
c.dat$period= c.dat$Time
c.dat$period[c.dat$period=="historic"]= "1961-1971"
c.dat$period[c.dat$period=="recent"]= "2001-2011"

#beta fits
co.colias= temp.co + geom_line(data=ps.all[ps.all$Population=="CO",], aes(x=Temp, y=Perf/5, color=period),size=1.1)
#divide by 5 to line up, use shade temps? 
ca.colias= temp.ca + geom_line(data=ps.all[ps.all$Population=="CA",], aes(x=Temp, y=Perf/5, color=period),size=1.1)

#add data
co.colias= co.colias + geom_point(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, y=mean/5, color=period),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5, color=period))+
  ylab("Feeding rate (g/g/h)")

ca.colias= ca.colias + geom_point(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, y=mean/5, color=period),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5, color=period))+
  ylab("Feeding rate (g/g/h)")

#reference temperatures
co.colias.ref= temp.co.ref + geom_line(data=ps.all[ps.all$Population=="CO",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)+
  geom_point(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, y=mean/5, shape=Time),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CO",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5))
ca.colias.ref= temp.ca.ref + geom_line(data=ps.all[ps.all$Population=="CA",], aes(x=Temp, y=Perf/5, lty=Time),size=1.1)+
  geom_point(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, y=mean/5, shape=Time),size=3)+
  geom_errorbar(data=c.dat[c.dat$Pop=="CA",], aes(x=Temp, ymin=mean/5-se/5, ymax=mean/5+se/5))

#==============================
#Colias feeding MARCH

#Fig 1: modern
#Fig 3: Historic, Sacramento and Montrose; 1961-1971, 2001-2011 
#analyze feeding rate distributions?

#plot temperature distributions
locations= c("Montrose","Sacramento", "LosBanos")

loc.k=2

#years for data
years=c(1961:1968,1971,2001:2010, 2011:2016,2018:2021) 

#combine data
for(yr.k in 1:length(years)){
  dat= read.csv(paste("./data/era5_micro_shade/",locations[loc.k],years[yr.k],".csv",sep="") )
  dat$year= years[yr.k]
  
  dat$period<- NA
  #vary periods
  if(years[yr.k]<1972) dat$period<- "1961-1971"
 # if(years[yr.k]>=2000 & years[yr.k]<=2010) dat$period <- "2000-2010"
  if(years[yr.k]>=2011 & years[yr.k]<=2021) dat$period <- "2011-2021"
  
  if(yr.k==1) dat.all=dat
  if(yr.k>1) dat.all=rbind(dat, dat.all)
}

#keep all data, no longer subset to sunlight
dat.day=dat.all

#Recode as Apr + May, July +Aug; 91:151, 182:243
dat.day$seas= NA
#April to May
dat.day$seas[dat.day$DOY %in% 60:120]= "MarApr"
# July to August
dat.day$seas[dat.day$DOY %in% 182:243]= "JulAug"
# #order
dat.day$seas= factor(dat.day$seas, levels=c("MarApr","JulAug") )
#drop other days
dat.day= dat.day[which(!is.na(dat.day$seas)),]

#combine doy and time
dat.day$d.hr= dat.day$DOY + dat.day$TIME/60

#drop middle period
dat.day.plot= dat.day[!is.na(dat.day$period),]
#make character factor
dat.day.plot$period= as.factor(as.character(dat.day.plot$period))

#make character factor
dat.day.plot$period= as.factor(as.character(dat.day.plot$period))

#plot density distributions
p1n= ggplot(dat.day.plot, aes(x=TALOC))+
  geom_density(alpha=0.3, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Density")+ #"Feeding rate (g/g/h)"
  xlab("Temperature (°C)" )+
  xlim(-5,50)+
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,0.075), expand=c(0,0)) +
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 
#D0cm, TALOC, TAREF

#plot reference temperatures
p1n.ref= ggplot(dat.day.plot, aes(x=TAREF))+
  geom_density(alpha=0.3, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Density")+
  xlab("Temperature at reference height (°C)" )+
  xlim(-5,50)+ 
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.6, 0.7),legend.background = element_rect(fill="transparent"))+
  scale_y_continuous(limits=c(0,0.075), expand=c(0,0)) +
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 

#plot together
temp.nielsen=p1n + geom_density(alpha=0.3, linetype="dashed", aes(x=TAREF, color=period))+
  theme(legend.position = c(0.5, 0.85),legend.background = element_rect(fill="transparent")) 

#--------------------
# Selection analysis

#a= height; b= mode; c= breadth; d= inverse breadth?; e= breadth

pop.times= c("CO_historic","CA_historic")

for(loc.k in 1:length(pop.times)){

#generate parameter combinations
tpc.beta= as.numeric(tpc.betas[which(tpc.betas$pop==pop.times[loc.k]),c(1:4,6)] ) 

param.grid= expand.grid(shift= seq(tpc.beta[1]-15,tpc.beta[1]+10,1),
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
perfs1= aggregate(perfs[,4:ncol(perfs)], list(perfs$year,perfs$seas), FUN=mean)
names(perfs1)=c("year","seas", 1:nrow(param.grid) )

#estimate topt for curve
param.grid$Topt=NA
for(topt.k in 1:nrow(param.grid)){
  ps= TPC_beta(1:50,shift=param.grid[topt.k,1], breadth=param.grid[topt.k,2], aran=0, 
               tolerance=tpc.beta[3], skew=tpc.beta[4])
  param.grid$Topt[topt.k]= c(1:50)[which.max(ps)]
}

#plot fitness curves through time
perfs.l= melt(perfs1[,1:ncol(perfs1)], id.vars = c("year","seas")) #change from perfs1[,1:(ncol(perfs)-1)]
names(perfs.l)[3:4]=c("temperature","performance")

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
fig.fitnesscurves=ggplot(perfs.l, aes(x=topt, y=performance, color=year, group=yrbr) )+geom_line(aes(lty=breadth),linewidth=1)+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+ xlim(-5,50)+
  xlab("Thermal optima (C)")+ylab("Feeding rate (g/g/h)") #+theme(legend.position = c(0.8, 0.3))

#fitness at fixed breadth
fig.fit.fb=ggplot(perfs.l[perfs.l$breadth==0.15,], aes(x=topt, y=performance, color=year, group=year) )+geom_line(linewidth=1)+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+ xlim(-5,50)+
  xlab("Thermal optima (C)")+ylab("Feeding rate (g/g/h)")+
  theme(legend.position = "right")+ #theme(legend.position = c(0.6, 0.8))+
  theme(strip.text.x = element_blank())

#fitness at fixed shift
fig.fit.fs=ggplot(perfs.l[round(perfs.l$shift,3)== 3.652,], aes(x=breadth, y=performance, color=year, group=year) )+geom_line()+
  facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+ xlim(-5,50)+
  xlab("Thermal optima (C)")+ylab("Mean feeding rate (g/g/h)")+
  theme(legend.position = "right") #theme(legend.position = c(0.8, 0.3))

#find thermal optima through years  
breadths= sort(unique(param.grid$breadth))

inds= which(param.grid$breadth==breadths[1])+2
perfs1$shift_opt_b1= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

inds= which(param.grid$breadth==breadths[2])+2
perfs1$shift_opt_b2= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

inds= which(param.grid$breadth==breadths[3])+2
perfs1$shift_opt_b3= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

#melt 
perfs.b= perfs1[,c("year","seas","shift_opt_b1","shift_opt_b2","shift_opt_b3")]
perfs.b.l= melt(perfs.b, id.vars = c("year","seas"))
names(perfs.b.l)[3:4]=c("breadth","opt_shift")
perfs.b.l$breadths= NA
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b1"]<- breadths[1]
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b2"]<- breadths[2]
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b3"]<- breadths[3]

perfs.b.l$seas_br= paste(perfs.b.l$seas, perfs.b.l$breadths,sep="_")
perfs.b.l$breadth= factor(perfs.b.l$breadths)

#beta= 0.15
fig.shift_opt= ggplot(perfs.b.l[which(perfs.b.l$breadth==0.15),], aes(x=year, y=opt_shift, color=seas, group= seas_br, lty=breadth))+geom_line(linewidth=1)+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  xlab("Year")+ylab("Optimal thermal optima (C)")+
  guides(lty="none")+
  theme(legend.position = "right")+ #theme(legend.position = c(0.2, 0.9))+ 
  scale_color_manual(values=c("#39568CFF", "#DCE319FF"))+
 # theme(legend.position = "bottom")+
  labs(colour = "season") 

#all betas
fig.shift_opt.all= ggplot(perfs.b.l, aes(x=year, y=opt_shift, color=seas, group= seas_br, lty=breadth))+geom_line(linewidth=1)+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  xlab("Year")+ylab("Optimal thermal optima (C)")+
  scale_color_manual(values=c("#39568CFF", "#DCE319FF"))+
  labs(colour = "season") 

if(loc.k==1) {fig.fitnesscurves.co= fig.fit.fb; fig.fit.all.co= fig.fitnesscurves}
if(loc.k==2) {fig.fitnesscurves.ca= fig.fit.fb; fig.fit.all.ca= fig.fitnesscurves}

if(loc.k==1) {fig.shift_opt.co= fig.shift_opt; fig.shift_opt.all.co= fig.shift_opt.all}
if(loc.k==2) {fig.shift_opt.ca= fig.shift_opt; fig.shift_opt.all.ca= fig.shift_opt.all}

} # end loop locations

#-------------------
#Feeding rate distributions and sums through time

pop.hist= c("CO_historic","CA_historic")
pop.rec= c("CO_recent","CA_recent")

for(loc.k in 1:2){
  
  #generate parameter combinations
  tpc.hist= as.numeric(tpc.betas[which(tpc.betas$pop==pop.hist[loc.k]),c(1:4,6)] ) 
  tpc.rec= as.numeric(tpc.betas[which(tpc.betas$pop==pop.rec[loc.k]),c(1:4,6)] ) 
 
  if(loc.k==1) dat.day= dat.day.co
  if(loc.k==2) dat.day= dat.day.ca
  
  #estimate performance
  dat.day$perf.histTPC= TPC_beta(dat.day$TALOC, shift=tpc.hist[1], breadth=tpc.hist[2], aran=0, tolerance=tpc.hist[3], skew=tpc.hist[4])/(2.5/tpc.hist[5])
  dat.day$perf.recTPC= TPC_beta(dat.day$TALOC, shift=tpc.rec[1], breadth=tpc.rec[2], aran=0, tolerance=tpc.rec[3], skew=tpc.rec[4])/(2.5/tpc.rec[5])
  
  #plot density distributions of performance
  dat.perf= dat.day[,c("period","seas","perf.histTPC","perf.recTPC")]
  #to long format
  dat.perf= gather(dat.perf, TPC, performance, perf.histTPC:perf.recTPC, factor_key=TRUE)
  
  perf.dens= ggplot(dat.perf, aes(x=performance))+
    geom_density(alpha=0.5, aes(fill=period, color=period, lty=TPC, linewidth=1))+
    facet_wrap(~seas)+
    xlab("Performance")+
    ylab("Density" )+
    theme_classic(base_size = 20)+theme(legend.position = c(0.5, 0.8))+
    ylim(0,0.75)+
    scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
    scale_fill_manual(values=c("#7AD151FF","#440154FF")  )  
  
  #performance means
  dat.day1= dat.day #na.omit(dat.day)
  perf.mean=cbind(aggregate(dat.day1$perf.histTPC, list(dat.day1$period, dat.day1$seas), FUN=mean),
  aggregate(dat.day1$perf.recTPC, list(dat.day1$period, dat.day1$seas), FUN=mean)[,3] )
  names(perf.mean)= c("period","season","histTPC","recTPC")
  
  #density distribution by performance
  dat.day1$Tround= round(dat.day1$TALOC)
  dat.day2= aggregate(dat.day1[,c("perf.histTPC","perf.recTPC")], list(dat.day1$Tround, dat.day1$period, dat.day1$seas), FUN=sum)
  names(dat.day2)= c("temp","period","seas","sumperf.hist","sumperf.rec")
  dat.day2$per_seas= paste(dat.day2$period, dat.day2$seas, sep="_")
  
  #normalize to count of data
  counts= aggregate(dat.day1$perf.histTPC, list(dat.day1$period, dat.day1$seas), FUN=function(x)length(x))
  names(counts)= c("period","seas","count")
  counts$per_seas= paste(counts$period, counts$seas, sep="_")
  #add counts to performance data
  dat.day2$counts= counts$count[match(dat.day2$per_seas, counts$per_seas)]
  #normalize
  dat.day2$sumperf.hist.norm= dat.day2$sumperf.hist/dat.day2$counts
  dat.day2$sumperf.rec.norm= dat.day2$sumperf.rec/dat.day2$counts
  
  #to long format
  names(dat.day2)[8:9]=c("historicTPC","recentTPC")
  dat.perf= gather(dat.day2[,c(1:3,8:9)], TPC, performance, historicTPC:recentTPC, factor_key=TRUE)
  
  feedt.hist= ggplot(dat.perf, aes(x=temp, y=performance, color=period, lty=TPC))+ #geom_line()+   
    geom_smooth(se=FALSE)+
    facet_wrap(~seas)+
    xlab("Temperature at reference height (°C)")+
    ylab("Sum of feeding rate (g/g/h)" )+
    #ylim(0,0.0015)+
    xlim(0,40)+
    ylim(0,0.75)
    theme_classic(base_size = 20)+theme(legend.position = c(0.5, 0.7))
  
  #----
  #performance over time for historic and recent TPCs
  perf.yr=cbind(aggregate(dat.day1$perf.histTPC, list(dat.day1$year, dat.day1$seas), FUN=mean),
                  aggregate(dat.day1$perf.recTPC, list(dat.day1$year, dat.day1$seas), FUN=mean)[,3] )
  names(perf.yr)= c("year","season","initial","recent")
  
  #to long format
  perf.yr.l= melt(perf.yr, id.vars = c("year","season"), variable.name = "tpc")
  names(perf.yr.l)[3]<-"TPC"
  
  tpc.yr= ggplot(perf.yr.l, aes(x=year, y=value, color=TPC))+ geom_line(linewidth=1)+   
    geom_smooth(method="lm",se=FALSE)+
    facet_wrap(~season)+
    xlab("Year")+
    ylab("Mean feeding rate (g/g/h)" )+
    theme_classic(base_size = 20)+theme(legend.position ="right")+ #c(0.7, 0.3)
    scale_color_manual(values=c("#7AD151FF","#440154FF"))
    #remove label
    #theme(strip.text.x = element_blank())
  
  #----
  #percent of time over thermal optima
  
 #find Topt
  hist.tpc<- tpc.betas[c(2,1),]
  rec.tpc<- tpc.betas[c(4,3),]
  
  ps= TPC_beta(1:50,shift=hist.tpc[loc.k,1], breadth=hist.tpc[loc.k,2], aran=0, 
               tolerance=hist.tpc[loc.k,3], skew=hist.tpc[loc.k,4])
  Topt.hist= c(1:50)[which.max(ps)]
  
  ps= TPC_beta(1:50,shift=rec.tpc[loc.k,1], breadth=rec.tpc[loc.k,2], aran=0, 
               tolerance=rec.tpc[loc.k,3], skew=rec.tpc[loc.k,4])
  Topt.rec= c(1:50)[which.max(ps)]
  
  tpc.betas$Ctmax<- tpc.betas$shift + tpc.betas$tolerance
  tpc.betas$t80<- tpc.betas$Ctmax - tpc.betas$tolerance*0.2
 # ctmax.hist<- tpc.betas$Ctmax[c(2,1)[loc.k]]-(tpc.betas$Ctmax[c(2,1)[loc.k]]-Topt.hist)/4
 # ctmax.rec<- tpc.betas$Ctmax[c(4,3)[loc.k]]-(tpc.betas$Ctmax[c(2,1)[loc.k]]-Topt.rec)/4
  ctmax.hist<- tpc.betas$t80[c(2,1)[loc.k]]
  ctmax.rec<- tpc.betas$t80[c(4,3)[loc.k]]
  
  #proportion above Topt
  #dat.day1$ext.hist<- ifelse(dat.day1$TALOC>Topt.hist, 1, 0)
  #dat.day1$ext.rec<- ifelse(dat.day1$TALOC>Topt.rec, 1, 0)
  
  dat.day1$ext.hist<- ifelse(dat.day1$TALOC>ctmax.hist, 1, 0)
  dat.day1$ext.rec<- ifelse(dat.day1$TALOC>ctmax.rec, 1, 0)
  
  #performance over time for historic and recent TPCs
  ext.yr=cbind(aggregate(dat.day1$ext.hist, list(dat.day1$year, dat.day1$seas), FUN=mean),
                aggregate(dat.day1$ext.rec, list(dat.day1$year, dat.day1$seas), FUN=mean)[,3] )
  names(ext.yr)= c("year","season","initial","recent")
  
  #proportion exceedences
  ext.yr[ext.yr$year %in% c(1972, 2012),]
  
  #to long format
  ext.yr.l= melt(ext.yr, id.vars = c("year","season"), variable.name = "tpc")
  names(ext.yr.l)[3]<-"TPC"
  
  ext.yr.plot= ggplot(ext.yr.l, aes(x=year, y=value, color=TPC))+ geom_line(linewidth=1)+   
    geom_smooth(method="lm",se=FALSE)+
    facet_wrap(~season)+
    xlab("Year")+
    ylab("Proportion extremes" )+
    theme(legend.position = "none")+
    theme_classic(base_size = 20)+theme(legend.position = "right")+ #c(0.2, 0.8)
    scale_color_manual(values=c("#7AD151FF","#440154FF"))+
    #remove label
    theme(strip.text.x = element_blank())
  
  #----
  #save
  if(loc.k==1) {feedfig.co= feedt.hist; perf.co=perf.mean; perf.dens.co= perf.dens; feed.tpc.co= tpc.yr; ext.yr.co= ext.yr.plot}
  if(loc.k==2) {feedfig.ca= feedt.hist; perf.ca=perf.mean; perf.dens.ca= perf.dens; feed.tpc.ca= tpc.yr; ext.yr.ca= ext.yr.plot}
  
} # end loop locations
  
  #combine performance means
  perf.co$location="CO"
  perf.ca$location="CA"
  perf.tpcs= rbind(perf.co, perf.ca)
  
  #to long format
  perf.mean= gather(perf.tpcs, TPC, performance, histTPC:recTPC, factor_key=TRUE)
  perf.mean$location= factor(perf.mean$location, levels=c("CO","CA"), ordered=TRUE)
  
  fig.meanperf= ggplot(perf.mean, aes(x=TPC, y=performance, color=period, shape=season))+
    facet_wrap(~location)+
    geom_point(size=3)+
    #scale_shape_manual(values=c(21,22))+
    ylab("mean performance")

 pdf("./figures/Fig_Colias_FeedingRateByTemp.pdf", height = 10, width = 10)
  feedfig.co / feedfig.ca
  dev.off()
  
  pdf("./figures/Fig_Colias_FeedingRateDensity.pdf", height = 10, width = 10)
  perf.dens.co / perf.dens.ca
  dev.off()
  
  pdf("./figures/MeanPerformance.pdf", height = 6, width = 10)
  fig.meanperf
  dev.off()
  
  pdf("./figures/Fig_Colias_PerfByYr.pdf", height = 10, width = 10)
  feed.tpc.co / feed.tpc.ca |
    ext.yr.co/ ext.yr.ca
  dev.off()
  
#-------------------

#PLOT

  design <- c("
    11333
    11333
    11444
    22444
    22555
    22555
  ")  
  
pdf("./figures/Fig2_Colias_CO_sun.pdf",height = 12, width = 20)
temps.co/
  co.colias/
  #fig.fitnesscurves.co/
  feed.tpc.co/
  ext.yr.co/
  fig.shift_opt.co +
  plot_layout(design = design)+
  plot_annotation(tag_levels = 'A')
dev.off()

pdf("./figures/Fig3_Colias_CA_sun.pdf",height = 12, width = 20)
temps.ca/
  ca.colias/
  #fig.fitnesscurves.ca/
  feed.tpc.ca/
  ext.yr.ca/
  fig.shift_opt.ca +
  plot_layout(design = design)+
  plot_annotation(tag_levels = 'A')
dev.off()

#supplementary plots
pdf("./figures/Fig_Colias_CO_supp_sun.pdf",height = 10, width = 8)
fig.fit.all.co/
  fig.shift_opt.all.co +
  plot_annotation(tag_levels = 'A')
dev.off()

pdf("./figures/Fig_Colias_CA_supp_sun.pdf",height = 10, width = 8)
fig.fit.all.ca/
  fig.shift_opt.all.ca +
  plot_annotation(tag_levels = 'A')
dev.off()

#-----
pdf("./figures/Fig5_photo_sun.pdf",height = 10, width = 8)
temp.nielsen/ fig.co.photo +
  plot_annotation(tag_levels = 'A')
dev.off()

#----
#save TPC parameters
tpc.betas$pop<- as.character(tpc.betas$pop)

#add Pieris
tpc.pr<- tpc.betas[1,]
tpc.pr[]<- NA
tpc.pr[1:4]<- c(-3.4132807, 0.1500000,  44.9612415, 0.7321506)
tpc.pr$pop<- "WA_historic"
tpc.pr$Ctmax<- tpc.pr$shift + tpc.pr$tolerance
tpc.betas.all<- rbind(tpc.betas, tpc.pr)

tpc.betas.all$Topt<-NA
tempv<- seq(from=1, to=50, by=0.1)

for(loc.k in 1:5){
ps= TPC_beta(tempv,shift=tpc.betas.all[loc.k,1], breadth=tpc.betas.all[loc.k,2], aran=0, 
             tolerance=tpc.betas.all[loc.k,3], skew=tpc.betas.all[loc.k,4])
tpc.betas.all$Topt[loc.k]<- tempv[which.max(ps)]
}


