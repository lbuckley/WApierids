## TO DO
# Add size into fitness estimates, P. rapae 
#Fig 5 in Jones et al. life time fecundity, y= 4.84w -273, w=weight(mg)

#CT ovipositing https://doi.org/10.1890/05-0647
#Pieris virginiensis
#Table 1, 0.36 host plants investigated per minute
#egg survival= 0.57
#larval survival= exp(0.586-0.153*day)/(1+exp(0.586-0.153*day))
larv.surv= function(day) exp(0.586-0.153*day)/(1+exp(0.586-0.153*day))

#Egg laying in Pieris rapae, https://www.jstor.org/stable/pdf/2401827.pdf
#temperature dependent survival rate
#egg production as a function of temperature accumulation during day
#Fig 6: proportion adults surviving to age x (degree days greater than 10C)= 1/(1+0.013*exp(0.28*x))
add.surv= function(dd10) 1/(1+0.013*exp(0.028*dd10)) #0.028 changed from 0.28 in paper to match plot

#https://www.jstor.org/stable/3956
#P. rapae vancouver
#eggs per visit: cabbage= 0.68, kale= 0.49, radish= 0.39
#use 0.52?

#https://www.jstor.org/stable/pdf/4217728.pdf
#Davies and Gilbert 1985
#P. rapae larval development, Ascot, UK
#eggs: 10C developmental threshold, 54 degree days to hatch
#larval: 157 DD_10 for larval development
#pupal development: 116 DD_9.3
# 18% egg survival

#https://www.jstor.org/stable/pdf/4536.pdf
#P. rapae, Vancouver
#size heritability
#survival data

#fit absolute growth rate
#https://doi.org/10.1111/j.0014-3820.2004.tb01732.x
#Fig 2, mg/hr for 4th instaar
temps= c(11,17,23,29,35,40)
growth= c(0.192,0.448,0.993,1.111,1.461,0.700)
#or estimate pupal weight as a function of temperature, Fig 2,Control of Fecundity II
growth_mgh= function(temp) weibull_1995(temp, a=1.415,topt=32.6,b=4.12*10^5,c=5.27*10^4)

#convert from wing traits to absorptivity

library(truncnorm)
library(dplyr)
library(TrenchR)
library(ggplot2)
library(MCMCglmm)
library(reshape2)
library(patchwork)
library(rTPC)

#LOAD PARAMETERS
#Demographic parameters
#Kingsolver 1983 Ecology 64
#OviRate=0.73 #Ovipositing rates: 0.73 eggs/min (Stanton 1980) 

#P. rapae vancouver, https://www.jstor.org/stable/3956
#eggs per visit: cabbage= 0.68, kale= 0.49, radish= 0.39, use mean 0.52 eggs per visit 
#CT ovipositing https://doi.org/10.1890/05-0647
#Pieris virginiensis, Table 1, 0.36 host plants investigated per minute
OviRate= 0.36*0.52 #eggs/minute
  
#MaxEggs=700; # Max egg production: 700 eggs (Tabashnik 1980)
PropFlight= 0.5; # Females spend 50% of available activity time for oviposition-related
# Watt 1979 Oecologia

#SurvDaily=0.6; # Daily loss rate for Colias p. eriphyle at Crested Butte, female values
#Hayes 1981, Data for Colias alexandra at RMBL
#SurvMat=0.014; #1.4# survival to maturity

#read/make species data
solar.abs= 0.65 # Solar absorptivity, proportion
SpecDat= as.data.frame(solar.abs)
SpecDat$d=0.26 # d- Thoractic diameter, cm
SpecDat$fur.thickness=1.46 # Fur thickness, mm
SpecDat= as.matrix(SpecDat)

#Fur thickness from Kingsolver (1983)
#Eriphyle Montrose 0.82mm, d=3.3mm, 1.7km
#Eriphyle Skyland 1.08, d=3.5mm, 2.7km
#Meadii Mesa sco 1.46mm, d=3.6mm, 3.6km
FT= 1.46 #c(0.01, 0.82, 1.46, 2)
#FUR THICKNESS: 0, eriphyle olathe, meadii mesa seco, 2

#------------
#Load envi data

locations= c("Corfu","Seattle")
loc.k=1

#years for data
if(loc.k==1) years=c(1989:2021) #1989:1993, 2017:2021
if(loc.k==2) years=c(2001:2021)

setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/')
#combine data
for(yr.k in 1:length(years)){
  dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
  dat$year= years[yr.k]
  if(years[yr.k]<2000)dat$period="1989-1999"
  if(years[yr.k]>=2000 & years[yr.k]<2011)dat$period="2000-2010"
  if(years[yr.k]>=2011)dat$period="2011-2021"
  
  if(yr.k==1) dat.all=dat
  if(yr.k>1) dat.all=rbind(dat, dat.all)
}

#subset to sunlight
dat.day= subset(dat.all, dat.all$SOLR>0)

#divide by periods, 91 to 274
#April to May
dat.day$seas="Apr & May"
#June to July
dat.day$seas[dat.day$DOY>151]="Jun & Jul"
#August to September
dat.day$seas[dat.day$DOY>212]="Aug & Sep"
#order
dat.day$seas= factor(dat.day$seas, levels=c("Apr & May","Jun & Jul","Aug & Sep") )

#combine doy and time
dat.day$d.hr= dat.day$DOY + dat.day$TIME/60

#subset columns
dat.sub= dat.day[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","ZEN","SOLR","VLOC","D0cm","year","period","seas")]
#TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
#TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
#VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)

#partition solar radiation [or extract from micro?], returns diffuse fraction
df=partition_solar_radiation("Erbs", kt=0.7) #0.8

#---------
#Load shade data

setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_shade/')
#combine data
for(yr.k in 1:length(years)){
  dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
  dat$year= years[yr.k]
  if(years[yr.k]<2000)dat$period="1989-1999"
  if(years[yr.k]>=2000 & years[yr.k]<2011)dat$period="2000-2010"
  if(years[yr.k]>=2011)dat$period="2011-2021"
  
  if(yr.k==1) dat.all.sh=dat
  if(yr.k>1) dat.all.sh=rbind(dat, dat.all)
}

#subset to sunlight
dat.day.sh= subset(dat.all.sh, dat.all.sh$SOLR>0)

#divide by periods, 91 to 274
#April to May
dat.day.sh$seas="Apr & May"
#June to July
dat.day.sh$seas[dat.day.sh$DOY>151]="Jun & Jul"
#August to September
dat.day.sh$seas[dat.day.sh$DOY>212]="Aug & Sep"
#order
dat.day.sh$seas= factor(dat.day.sh$seas, levels=c("Apr & May","Jun & Jul","Aug & Sep") )

#combine doy and time
dat.day.sh$d.hr= dat.day.sh$DOY + (dat.day.sh$TIME/60-7)/24

#subset columns
dat.sub.sh= dat.day.sh[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","ZEN","SOLR","VLOC","D0cm","year","period","seas")]
#TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
#TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
#VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)

#match
match1= match(dat.sub$dates, dat.sub.sh$dates)

#compare data
plot(dat.sub.sh$TALOC[match1],dat.sub$TALOC)
plot(dat.sub.sh$D0cm[match1],dat.sub$D0cm)

#add shade data to sun data
dat.sub$TALOC_sh= dat.sub.sh$TALOC[match1]
dat.sub$D0cm_sh= dat.sub.sh$D0cm[match1]

#plot time series
#melt temp data
datm= melt(dat.sub[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","D0cm","year","period")], id=c("dates","DOY","TIME","d.hr","year","period") )

#plot
ggplot(data=datm, aes(x=d.hr, y = value, color=variable))+ geom_line(alpha=0.2)+
  theme_bw()+scale_color_viridis_d()

#----------
#butterfly temperature
Temat= cbind(dat.sub$TALOC, dat.sub$TALOC_sh,dat.sub$D0cm, dat.sub$D0cm_sh, 
             dat.sub$VLOC, dat.sub$SOLR*(1-df), dat.sub$SOLR*(df), dat.sub$ZEN)

dat.sub$Tb= apply(Temat, MARGIN=1, FUN=Tb_butterfly.mat,
                  D = 0.26, delta = 1.46, HB = 0.55, PV = 0.55, 
                  r_g = 0.3, wing_angle=42, shade = TRUE) 

#fix high
dat.sub$Tb[which(dat.sub$Tb>80)]<-NA

#plot Tb distributions
ggplot(data=dat.sub, aes(x=d.hr, y = Tb))+ geom_line(alpha=0.4)+
  theme_bw()+scale_color_viridis_d()
#plot density distributions
p1= ggplot(dat.sub, aes(x=Tb))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Body Temperature (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.65, 0.85))+
  xlim(0,90)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccTb.pdf", height = 7, width = 10)
p1
dev.off()

#air temp
p2= ggplot(dat.sub, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Environmental Temperature (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.7, 0.85))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccTaloc.pdf", height = 7, width = 10)
p2
dev.off()

#plot diurnal variation
#seems to work but still not sure about time
dat.sub$hr= dat.sub$TIME/60+12+7
dat.sub$hr[dat.sub$hr<0]= 24+dat.sub$hr[dat.sub$hr<0]
dat.sub$hr[dat.sub$hr>24]= dat.sub$hr[dat.sub$hr>24]-24

ggplot(dat.sub, aes(x=hr, y=Tb))+
  geom_point()+  
  facet_wrap(~seas)

#---------
### SET UP DATA STRUCTURES
#Make array to store data
#Lambda= survival*fecundity
#Matrix of lambdas
#dims: stats, year, Lambda 
Lambda<-array(NA, dim=c(length(years),length(seq(0.4,0.7,0.05)),length(seq(0.4,0.7,0.05)),5,5)) 
#2nd and 3rd dimensions are PV and HB
#Last dimension is Lambda, FAT,Egg Viability, growth
dimnames(Lambda)[[1]]<-years
dimnames(Lambda)[[4]]<-c("gen1","gen2","gen3","gen4","gen5")

#Matrix for pupual temps
pup.temps<-array(NA, dim=c(12, length(years),5)) #5 generations 
#Add names
dimnames(pup.temps)[[1]]= c("stat","yr","gen","Jlarv", "Jpup","Jadult","Tlarv","Tpup","Tad","Tlarv_fixed","Tpup_fixed","Tad_fixed")

#LOOP YEARS
years=c(1989:2021)

for(yr.k in 1:length(years) ){
  print(yr.k)
  
  #subset to year
  dat.yr= subset(dat.sub, dat.sub$year==years[yr.k])
  #need all days in year
  
  #============================================================================
  #CALCULATE DEVELOPMENT TIMING AND TEMPS
 
  #https://www.jstor.org/stable/pdf/4217728.pdf
  #Davies and Gilbert 1985
  #P. rapae larval development, Ascot, UK
  #eggs: 10C developmental threshold, 54 degree days to hatch
  #larval: 157 DD_10 for larval development
  #pupal development: 116 DD_9.3
  
  # ESTIMATE DEVELOPMENTAL TIMING
  #https://www.jstor.org/stable/4913
  
  DevZeros= 10 #egg to larval
  DevZeros_pup= 9.3 #pupal
  
  GddReqs_egg= c(157) #base 10C, https://www.jstor.org/stable/pdf/4217728.pdf
  GddReqs_larv= c(157) #base 10C, https://www.jstor.org/stable/pdf/4217728.pdf
  GddReqs_pup= c(157)  #base 9.3C, https://www.jstor.org/stable/pdf/4217728.pdf
  
  dat.yr$dd= dat.yr$TALOC -DevZeros
  dat.yr$dd[dat.yr$dd<0]= 0 #set negative to zero
  dat.yr$dd= dat.yr$dd/24 #divide by 24 to switch from degree days to hours
  dat.yr$dd.cs= cumsum(dat.yr$dd)
  
  dat.yr$dd_pup= dat.yr$TALOC -DevZeros_pup
  dat.yr$dd_pup[dat.yr$dd_pup<0]= 0 #set negative to zero
  dat.yr$dd_pup= dat.yr$dd_pup/24 #divide by 24 to switch from degree days to hours
  dat.yr$dd.cs_pup= cumsum(dat.yr$dd_pup)
  
  #estimate pupal timing
  #compare to empirical data
  setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/')
  pr= read.csv("PrapaeUW.Seln2.1999.Combineddata.OPUS2021.csv")
  
  #mean time to pupation and eclosion
  colMeans(pr[,c("Time.to.Pupation","Time.to.Eclosion")], na.rm=TRUE)
  #about 7 day pupation
  
  #-----------------------------------
  #Calculate development timing
  
    for(gen.k in 1:5){
      
      #Assume 7 days from eclosion to eggs laid ##UPDATE?
      #TO LARVAE
      Jegg.start= ifelse(gen.k>1, Jadult+7, dat.yr$DOY[which.max(dat.yr$dd>0)] )  
      if(Jegg.start>max(dat.yr$DOY)) Jegg.start=max(dat.yr$DOY)
      inds= which(dat.yr$DOY==Jegg.start)[1]:nrow(dat.yr)
      Jlarv= dat.yr$DOY[inds[which.max(cumsum(dat.yr$dd[inds])> GddReqs_egg)]] 
      if(Jlarv>max(dat.yr$DOY) | length(Jlarv)==0) Jlarv=max(dat.yr$DOY)
      
      ##TO PUPATION
      inds= which(dat.yr$DOY==Jlarv)[1]:nrow(dat.yr)
      Jpup= dat.yr$DOY[inds[which.max(cumsum(dat.yr$dd[inds])> GddReqs_larv)]] 
       if(Jpup>max(dat.yr$DOY) | length(Jpup)==0) Jpup=max(dat.yr$DOY) 
      
      #PUPATION 
      # ~7 day pupation based on CUH field data
      #Change to DD
      inds= which(dat.yr$DOY==Jpup)[1]:nrow(dat.yr)
      Jadult= dat.yr$DOY[inds[which.max(cumsum(dat.yr$dd_pup[inds])> GddReqs_pup)]] 
      if(Jadult>max(dat.yr$DOY) | length(Jadult)==0) Jadult=max(dat.yr$DOY)
      
      #----------------------
      #Calculate temps
      Tlarv= mean( dat.yr[dat.yr$DOY %in% Jlarv:(Jpup-1),"TALOC"] , na.rm=TRUE)
      Tpup= mean( dat.yr[dat.yr$DOY %in% Jpup:(Jadult-1),"TALOC"] , na.rm=TRUE)
      Tad= mean( dat.yr[dat.yr$DOY %in% Jadult:(Jadult+5),"TALOC"] , na.rm=TRUE)
      
      #Write data in array
      pup.temps[3:9,yr.k, gen.k]=c(gen.k,Jlarv,Jpup,Jadult,Tlarv,Tpup,Tad)
      
    } #end loop generation
} #end loop years
 
  #plot
pup.temps.l= rbind(t(pup.temps[,,1]),t(pup.temps[,,2]),t(pup.temps[,,3]),t(pup.temps[,,4]),t(pup.temps[,,5]) )
pup.temps.l= as.data.frame(pup.temps.l)
pup.temps.l$yr= rep(years,5)
pup.temps.l$gen= c(rep(1,33),rep(2,33),rep(3,33),rep(4,33),rep(5,33) )
  
pup.temps.l= melt(pup.temps.l[,2:9], id=c("yr", "gen")) 

#plot
plot.pocc.times= ggplot(pup.temps.l[which(pup.temps.l$variable %in% c("Jlarv","Jpup","Jadult") ),], aes(x=yr, y=value, col=gen, group=gen))+geom_line()+
  facet_wrap(~variable)+
  theme_classic(base_size = 20)+
  ylab("day of year")+
  scale_color_viridis_c("generation")
plot.pocc.temps= ggplot(pup.temps.l[which(pup.temps.l$variable %in% c("Tlarv","Tpup","Tad") ),], aes(x=yr, y=value, col=gen, group=gen))+geom_line()+
  facet_wrap(~variable)+
  theme_classic(base_size = 20)+
  ylab("temperature (C)")+
  scale_color_viridis_c("generation")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccPupTemps.pdf", height = 7, width = 10)
plot.pocc.times / plot.pocc.temps
dev.off()

  #=======================================
  #DEMOGRAPHY
  
#PERFORMANCE FUNCTIONS
#function for flight probability
fl.ph<- function(x) 1 * exp(-0.5*(abs(x-33.5)/5)^3.5) #Colias curve
#P. protodice, losely related to P. occidentalis, has similar curve: https://www.jstor.org/stable/pdf/4217668.pdf

#function for egg viability
egg.viab<-function(x) ifelse(x<40, egg.viab<-1, egg.viab<- exp(-(x-40)/35.32)) #colias

#define geometric mean
geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

  #-------------------------------------------
  #DEMOGRAPHY

#absorptivity
abs1=seq(0.4,0.7,0.05)

for(yr.k in 1:length(years) ){
  
    dat.yr= subset(dat.sub, dat.sub$year==years[yr.k])
  
    Temat= cbind(dat.yr$TALOC, dat.yr$TALOC_sh,dat.yr$D0cm, dat.yr$D0cm_sh, 
                 dat.yr$VLOC, dat.yr$SOLR*(1-df), dat.yr$SOLR*(df), dat.yr$ZEN)
    
  for(gen.k in 1:5 ){ #loop generation
    
      for(abs.k.pv in 1:length(abs1) ){ #loop PV absorptivity
        for(abs.k.hb in 1:length(abs1) ){ #loop HB absorptivity
          
        #ADD PLASTICITY
        #https://doi.org/10.1093/icb/38.3.545
        #Basal dorsal HW melanism
        #use females, check units of basal dorsal HW melanism
        #Adapt to HB, PV
        
        #hours light
        light.h= c(16,10)
        abs.m= c(.66,.67); abs.m= c(.717,.735)        
        mod1= lm(abs.m~light.h)
        abs.plast= function(light.hrs, abs) abs -0.003*light.hrs
        
        light.hrs=daylength(lat=46.8151, doy=pup.temps["Jpup", yr.k, gen.k])
        
        #absorptivity #UPDATE for pv, hb
        abs.pv= abs.plast(light.hrs, abs1[abs.k.pv])
        abs.hb= abs.plast(light.hrs, abs1[abs.k.hb])
        
        #calculate Tb based on absorptivity
        dat.yr$Tb.a= apply(Temat, MARGIN=1, FUN=Tb_butterfly.mat,
                          D = 0.26, delta = 1.46, HB = abs.hb, PV = abs.pv, 
                          r_g = 0.3, wing_angle=42, shade = TRUE) 
        
        #Flight probability
        dat.yr$fl.p= sapply(dat.yr$Tb.a, FUN=fl.ph)
        
        #Egg viability
        dat.yr$egg.v= sapply(dat.yr$Tb.a, FUN=egg.viab)
      
        #absolute growth rate
        #https://doi.org/10.1111/j.0014-3820.2004.tb01732.x
        #Fig 2, mg/hr for 4th instaar
        dat.yr$larv.growth= sapply(dat.yr$TALOC, FUN=growth_mgh)
        #aggregate larval growth as mean across hours and multiply by 24 to reach year
        
        #or estimate pupal weight as a function of temperature, Fig 2,Control of Fecundity II
        
       # daily values
        dat.day<- dat.yr %>%
          group_by(DOY) %>%
          summarise(FAT = sum(fl.p), EggViab= geo_mean(egg.v), Tb.mean= mean(Tb.a), Gt.l = mean(larv.growth)*24,
                    .groups="keep")
        
        #flight day
        Jfl= pup.temps["Jadult", yr.k, gen.k]

        ##CALCULATE EGG VIABILITY OVER 5 DAY PERIOD (GEOMETRIC MEAN ACROSS HOURS)
        #sample flight day from truncated normal distribution
        Nind=1000 #changed from 100
        f.low= max(Jfl-7,min(dat.sub$DOY)+2)
        f.up= min(Jfl+7,max(dat.sub$DOY)-2)
        
        flightday= round(rtruncnorm(Nind, a=f.low, b=f.up, mean = Jfl, sd = 2) )
        
        f.ind= match(flightday, dat.day$DOY)
        #if NA day, use mean
        f.ind[is.na(f.ind)]<-match(Jfl, dat.day$DOY)
        
        #calculate geometric mean of egg viability within flight period
        ev.ind=sapply(f.ind, function(x)  geo_mean(dat.day$EggViab[(x-2):(x+2)]) )
        #AVERAGE FAT OVER DAYS
        FAT.ind= sapply(f.ind, function(x)  mean(dat.day$FAT[(x-2):(x+2)], na.rm=TRUE) )
        #AVERAGE TEMP
        T.ind= sapply(f.ind, function(x)  mean(dat.day$Tb.mean[(x-2):(x+2)], na.rm=TRUE) )
        
        #LARVAL GROWTH
        ind.st= match(pup.temps["Jlarv",yr.k, gen.k], dat.day$DOY)
        ind.fin=match(pup.temps["Jpup",yr.k, gen.k], dat.day$DOY)
        
        Gr.larv= sum(dat.day$Gt.l[ind.st:(ind.fin-1)], na.rm=TRUE)
        
        Eggs.ind= 60*PropFlight*OviRate*FAT.ind * ev.ind #account for Egg viability
        Eggs.ind_noViab= 60*PropFlight*OviRate*FAT.ind
        
        #Means across individuals
        Eggs= mean(Eggs.ind)
        Eggs_noViab= mean(Eggs.ind_noViab)
        
        #Set max eggs as a function of ppupal size
        #https://doi.org/10.1071/ZO9820223
        #Fig 5 in Jones et al. life time fecundity, y= 4.84w -273, w=weight(mg)
        #bound pupal weight 110-210 mg
        if(Gr.larv<110) Gr.larv=110
        if(Gr.larv>210) Gr.larv=210
        MaxEggs= 4.84*Gr.larv-273
        
        #Survival to Maturity
        #egg surivival 0.18, Davies and Gilbert 1985, https://www.jstor.org/stable/pdf/4217728.pdf
        #pupal survival?
        days.larv= pup.temps["Jpup", yr.k, gen.k]-pup.temps["Jlarv", yr.k, gen.k]+1
        #P. virginiensis larval survival, https://doi.org/10.1890/05-0647
        SurvMat= 0.18*larv.surv(days.larv)
        
        if(!is.nan(Eggs)){
          MaxDay=5
          Lambda1=0
          for(day in 1:MaxDay){
            Eggs1= min(Eggs, MaxEggs-Eggs*(day-1))
            if(Eggs1<0) Eggs1=0
            
            #Pieris virginiensis, https://doi.org/10.1890/05-0647
            SurvDaily= mean(add.surv((T.ind-10)*day))
            
            Lambda1= Lambda1+ SurvMat * SurvDaily *Eggs1;                        
          }#end loop days
          
          Lambda[yr.k, abs.k.pv, abs.k.hb, gen.k, ]= c(Lambda1, mean(FAT.ind), mean(ev.ind), mean(T.ind, na.rm=T), Gr.larv )
          
        } #Check Eggs
        
      } #end loop pv absorptivity
        
      } #end loop hb absorptivity
     
  } #end loop generation
  
} #end loop years

#SAVE OBJECT
filename= paste("lambda1_",projs[proj.k],".rds",sep="")
saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#melt lambda array
Lambda.l= melt(Lambda) 
colnames(Lambda.l)=c("year","abs.pv","abs.hb","gen","metric","value")
Lambda.l$abs.pv= abs1[Lambda.l$abs.pv]
Lambda.l$abs.hb= abs1[Lambda.l$abs.hb]
Lambda.l$metric= c("fitness","FAT (h)","egg viab (%)","Tadult (C)","larval growth")[Lambda.l$metric]
Lambda.l$metric= factor(Lambda.l$metric, levels=c("fitness","FAT (h)","egg viab (%)","Tadult (C)","larval growth"))

#replace gen name
Lambda.l$gen= gsub("gen","generation ",Lambda.l$gen)

#plot lambdas
Lambda.l$gyr= paste(Lambda.l$gen, Lambda.l$year, sep="")

fig.OccFit.pv= ggplot(Lambda.l[which(Lambda.l$abs.hb==0.55 & Lambda.l$metric %in% c("fitness","FAT (h)","egg viab (%)")),], 
                   aes(x=abs.pv, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()+
  theme_classic()+xlab("PV absorptivity (%)")+ylab("")

fig.OccFit.hb= ggplot(Lambda.l[which(Lambda.l$abs.pv==0.55 & Lambda.l$metric %in% c("fitness","FAT (h)","egg viab (%)")),], 
                      aes(x=abs.hb, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()+
  theme_classic()+xlab("HB absorptivity (%)")+ylab("")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccFitCurv.pv.pdf", height = 7, width = 7)
fig.OccFit.pv
dev.off()

pdf("Fig_PoccFitCurv.hb.pdf", height = 7, width = 7)
fig.OccFit.hb
dev.off()

#just plot lambdas
fig.lambdas.pv= ggplot(Lambda.l[which(Lambda.l$abs.hb==0.55 & Lambda.l$metric %in% c("fitness")),], 
                      aes(x=abs.pv, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()+
  theme_classic()+xlab("PV absorptivity (%)")+ylab("fitness")

fig.lambdas.hb= ggplot(Lambda.l[which(Lambda.l$abs.pv==0.55 & Lambda.l$metric %in% c("fitness")),], 
                      aes(x=abs.hb, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()+
  theme_classic()+xlab("HB absorptivity (%)")+ylab("fitness")


pdf("Fig_PoccLambdas.pdf", height = 7, width = 7)
fig.lambdas.hb + fig.lambdas.pv + plot_layout(ncol = 1)
dev.off()

#-------------
#test Tb function
dat.sub1= dat.sub[dat.sub$DOY %in% c(100:102) & dat.sub$year==2021,]
Ts= matrix(NA, nrow=nrow(dat.sub1), ncol=length(abs1))

Temat= cbind(dat.sub1$TALOC, dat.sub1$TALOC_sh, dat.sub1$D0cm, dat.sub1$D0cm_sh, 
             dat.sub1$VLOC, dat.sub1$SOLR*(1-df), dat.sub1$SOLR*(df),dat.sub1$ZEN)

for(abs.k in 1:length(abs1) ){
  
  Ts[,abs.k]= apply(Temat, MARGIN=1, FUN=Tb_butterfly.mat,
                    D = 0.26, delta = 1.46, HB = abs1[abs.k], PV = abs1[abs.k], 
                    r_g = 0.3, wing_angle=42, shade = TRUE) 
  
  }

Ts= as.data.frame(Ts)
colnames(Ts)= abs1

Ts$DOY= dat.sub1$DOY
Ts$hr= dat.sub1$TIME/60 -12 #convert from UTC
Ts$hr[Ts$hr<0]= Ts$hr[Ts$hr<0]+24
Ts$TALOC= dat.sub1$TALOC
Ts= subset(Ts, Ts$DOY==101)
Ts= Ts[,c(1:7,9:10)]

#plot
datm= melt(Ts, id=c("hr"), variable="abs" )
ggplot(datm, aes(x=hr, y=value, color=abs, group=abs))+geom_line()

#========================================

#EVOLUTIONARY MODEL
N.ind=1000
a= seq(0.4,0.7,0.05)
#for finding a with max fitness
a.fit= as.data.frame(seq(0.4,0.7,0.01))
names(a.fit)="a"

##Read lambdas and pupal temps
##Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
#Lambda <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
#pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )

#Find years with calculations
counts= rowSums(is.na(pup.temps[,,1]))

inds=1:150
ngens=5

#==============================================
#Calculate optimal absorptivity

abs.opt= array(NA, dim=c(length(years), 5,2)) #last dimension is pv, hb  

lambda.fit=function(x, a=a) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients
lambda.max=function(x, a=a) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] 

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]
    
    #Extract temperatures
    Tp= pup.temps["Tpup",yr.k, gen.k]
    
    #--------------------------
    if(!all(is.na(Lambda.yr.gen[,,1]))){ #check has data
      
      inds= which(Lambda.yr.gen[,,1] == max(Lambda.yr.gen[,,1]), arr.ind = TRUE)
      abs.opt[yr.k,gen.k,1]= abs1[inds[1]] 
      abs.opt[yr.k,gen.k,2]= abs1[inds[2]]
      
      ##just hb
      #inds= which(Lambda.yr.gen[4,,1] == max(Lambda.yr.gen[4,,1]), arr.ind = TRUE)
      #abs.opt[yr.k,gen.k,1]= abs1[4] 
      #abs.opt[yr.k,gen.k,2]= abs1[inds]
      
    } #end check data
  } #end gen loop
} #end year loop

#save optimal Alphas
#saveRDS(abs.opt, paste("abs.opt_",projs[proj.k],".rds", sep=""))
#abs.opt <- readRDS( paste("abs.opt_",projs[proj.k],".rds", sep="") )

#plot optima alphas
abs.opt.l=as.data.frame(abs.opt[,,1])
abs.opt.l$year=years
abs.opt.l.pv=melt(abs.opt.l,id.vars = c("year"))
abs.opt.l.pv$trait="pv"

abs.opt.l=as.data.frame(abs.opt[,,2])
abs.opt.l$year=years
abs.opt.l.hb=melt(abs.opt.l,id.vars = c("year"))
abs.opt.l.hb$trait="hb"

abs.opt.l= rbind(abs.opt.l.pv, abs.opt.l.hb)

ggplot(abs.opt.l, aes(x=year, y=value, color=variable, lty=trait))+geom_line()+
  scale_color_viridis_d()

#***************************************
#compute initial AbsMean 
int_elev = 0.4226; slope_elev = 0.06517
Tmid = 20; slope_plast = -0.0083  #if Tmid=22.5, -0.006667;

h2.hb= (0.54+0.61)/2 #2 experiments #Kingsolver and Wiernasz 1991 Evolution
h2.pv= (0.22+0.23)/2 #2 experiments
cor=0.49

abs.sd= 0.062
rn.sd= 0.0083

#initialize with optimum value across for ten years, across generations
abs.init.pv <- mean(abs.opt[1:10,,1], na.rm=TRUE)
abs.init.hb <- mean(abs.opt[1:10,,2], na.rm=TRUE)

#-----------------------
#Save values
abs.mean= array(NA, dim=c(length(years), 5, 5,5,2))  #dims: yr.k, gen.k, scen.k:no plast, plast, only plast, metrics: abssample, absmid, rn, Babsmid, Brn; pv, hb
abs.mean[1,,,2,1]= abs.init.pv
abs.mean[1,,,3,1]= slope_plast
abs.mean[1,,,2,2]= abs.init.hb
abs.mean[1,,,3,2]= slope_plast
dimnames(abs.mean)[[4]]= c("abssample", "absmid", "rn", "Babsmid", "Brn") 

lambda.mean= array(NA, dim=c(length(years),5, 5,2)) #dims: yr.k, gen.k, scen.k:no plast, plast, only plast; pv, hb

#BetaRN= rep(NA, nrow(pts.sel))
#-------------------------------
scen.mat= rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,0),c(1,1,1) )
colnames(scen.mat)= c("plast","evol","evolRN"  )

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    BetaAbsmid.pv=NA
    BetaAbsmid.hb=NA
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k,]
    
    #determine those completing generations
    comp.gen= which(pup.temps["Jadult",yr.k,gen.k]<243)
    nocomp.gen= which(pup.temps["Jadult",yr.k,gen.k]==243)
    #set those not completing generations to NA
    if(length(nocomp.gen)>0) Lambda.yr.gen[nocomp.gen,nocomp.gen,]=NA
    
      #Extract temperatures
      Tp= pup.temps["Tpup",yr.k, gen.k]
      
      #Fitness models
      fit.pv= lm(Lambda.yr.gen[,4,1]~a+I(a^2))
      fit.hb= lm(Lambda.yr.gen[4,,1]~a+I(a^2))
      
      #Estimate fitness functions across cells
      fit.mod.pv= fit.pv$coefficients
      fit.mod.hb= fit.hb$coefficients
      
      #-------------------------
      # LOOP PLASTICITY SCENARIOS
      for(scen.k in 1:5){ #plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1
        
        if(scen.mat[scen.k,1]==1) rn.mean1= c(slope_plast,slope_plast)
        if(scen.mat[scen.k,1]==0) rn.mean1= c(0,0)
        if(scen.k==5 & gen.k==1) rn.mean1= abs.mean[yr.k,gen.k,scen.k,"rn",]
        if(scen.k==5 & gen.k>1) rn.mean1= abs.mean[yr.k,gen.k-1,scen.k,"rn",]
        
        if(gen.k==1) abs.mean1= abs.mean[yr.k,gen.k,scen.k,"absmid",]
        if(gen.k>1) abs.mean1= abs.mean[yr.k,gen.k-1,scen.k,"absmid",]
        
        #change NA values to negative values 
        abs.na.inds= abs.mean1[which( is.na(abs.mean1))]
        rn.na.inds= rn.mean1[which( is.na(rn.mean1))]
        
        abs.mean1[which( is.na(abs.mean1))]= -10 
        rn.mean1[which( is.na(rn.mean1))]= -1000
        
        #PV
        fit.mod=fit.mod.pv
        trait.k=1
        
        #Choose random sample of abs and rn values from current distribution (truncated normal) 
        abs.sample.pv= sapply(abs.mean1[trait.k], function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
        rn.sample.pv= sapply(rn.mean1[trait.k], function(x) rtnorm(N.ind, mean = x, sd = rn.sd, lower=-1, upper=1) )
        if(scen.mat[scen.k,1]==0) rn.sample.pv[]=0
        
        #Add plasticity across sites and sample
        abs.plast.pv <- abs.sample.pv + rn.sample.pv*(Tp-Tmid)
        #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
        
        ##calculate fitness 
        #use fitness function to predict Lambda for each individual
        #extract coefficients and calculate across abs samples
        fit.sample.pv=  sapply(abs.plast.pv[,1], function(x) fit.mod[1]+x*fit.mod[2]+x^2*fit.mod[3] )
       
        #standardize to relative fitness and centered on trait mean
        fit.mean.pv= mean(fit.sample.pv)
        lambda.mean[yr.k,gen.k,scen.k,trait.k]=fit.mean.pv
        rel.fit.pv= fit.sample.pv/fit.mean.pv
        
        absmid.dif.pv= t( apply(abs.sample.pv,1,'-',abs.mean1[trait.k]) )
        rn.dif.pv= t( apply(rn.sample.pv,1,'-',rn.mean1[trait.k]) )
        #---------
        #HB
        fit.mod=fit.mod.hb
        trait.k=2
        
        #Choose random sample of abs and rn values from current distribution (truncated normal) 
        abs.sample.hb= sapply(abs.mean1[trait.k], function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
        rn.sample.hb= sapply(rn.mean1[trait.k], function(x) rtnorm(N.ind, mean = x, sd = rn.sd, lower=-1, upper=1) )
        if(scen.mat[scen.k,1]==0) rn.sample.hb[]=0
        
        #Add plasticity across sites and sample
        abs.plast.hb <- abs.sample.hb + rn.sample.hb*(Tp-Tmid)
        #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
        
        ##calculate fitness 
        #use fitness function to predict Lambda for each individual
        #extract coefficients and calculate across abs samples
        fit.sample.hb=  sapply(abs.plast.hb[,1], function(x) fit.mod[1]+x*fit.mod[2]+x^2*fit.mod[3] )
        
        #standardize to relative fitness and centered on trait mean
        fit.mean.hb= mean(fit.sample.hb)
        lambda.mean[yr.k,gen.k,scen.k,trait.k]=fit.mean.hb
        rel.fit.hb= fit.sample.hb/fit.mean.hb
        
        absmid.dif.hb= t( apply(abs.sample.hb,1,'-',abs.mean1[trait.k]) )
        rn.dif.hb= t( apply(rn.sample.hb,1,'-',rn.mean1[trait.k]) )
        
        #---------
        
        R2selnAbsmid.pv<- 0 #No response to selection if no evolution
        R2selnRN.pv<- 0 
        R2selnAbsmid.hb<- 0 #No response to selection if no evolution
        R2selnRN.hb<- 0 
        #------------
        if(scen.k<5 & scen.mat[scen.k,2]==1){
          ##selection analysis
          #pv
          rel.fit.pv= as.vector(rel.fit.pv)
          absmid.dif.pv= as.vector(absmid.dif.pv)
          sel.fit.pv= lm(rel.fit.pv~absmid.dif.pv +I(absmid.dif.pv^2))$coefficients
          sel.mod.pv= lm(rel.fit.pv~absmid.dif.pv +I(absmid.dif.pv^2))
          BetaAbsmid.pv <-sel.fit.pv[2]
          
          #hb
          rel.fit.hb= as.vector(rel.fit.hb)
          absmid.dif.hb= as.vector(absmid.dif.hb)
          sel.fit.hb= lm(rel.fit.hb~absmid.dif.hb +I(absmid.dif.hb^2))$coefficients
          sel.mod.hb= lm(rel.fit.hb~absmid.dif.hb +I(absmid.dif.hb^2))
          BetaAbsmid.hb <-sel.fit.hb[2]
          
          #Response to selection
          R2selnAbsmid.pv <- h2.pv*(abs.sd^2)*BetaAbsmid.pv + cor*(abs.sd^2)*BetaAbsmid.hb
          R2selnAbsmid.hb <- h2.hb*(abs.sd^2)*BetaAbsmid.hb + cor*(abs.sd^2)*BetaAbsmid.pv
          
        } #end scen.k<5
        #------------
        if(scen.k==5){    
          ##selection analysis
          #pv
          rel.fit.pv= as.vector(rel.fit.pv)
          absmid.dif.pv= as.vector(absmid.dif.pv)
          sel.fit.pv= lm(rel.fit.pv~absmid.dif.pv +I(absmid.dif.pv^2))$coefficients
          sel.mod.pv= lm(rel.fit.pv~absmid.dif.pv +I(absmid.dif.pv^2))
          BetaAbsmid.pv <-sel.fit.pv[2]
          BetaRN.pv <- sel.fit.pv[3] 
          
          #hb
          rel.fit.hb= as.vector(rel.fit.hb)
          absmid.dif.hb= as.vector(absmid.dif.hb)
          sel.fit.hb= lm(rel.fit.hb~absmid.dif.hb +I(absmid.dif.hb^2))$coefficients
          sel.mod.hb= lm(rel.fit.hb~absmid.dif.hb +I(absmid.dif.hb^2))
          BetaAbsmid.hb <-sel.fit.hb[2]
          BetaRN.hb <- sel.fit.hb[3]
          
          #Response to selection
          R2selnAbsmid.pv <- h2.pv*(abs.sd^2)*BetaAbsmid.pv + cor*(abs.sd^2)*BetaAbsmid.hb
          R2selnAbsmid.hb <- h2.hb*(abs.sd^2)*BetaAbsmid.hb + cor*(abs.sd^2)*BetaAbsmid.pv
          
          R2selnRN.pv <- h2.pv*(rn.sd^2)*BetaRN.pv + cor*(rn.sd^2)*BetaRN.hb
          R2selnRN.hb <- h2.hb*(rn.sd^2)*BetaRN.hb + cor*(rn.sd^2)*BetaRN.pv
          
        } #end scen.k==5
        #-------------
        
        #Response to selection
        if(gen.k<5) {
          
          #PV
          trait.k=1
          abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]= abs.mean[yr.k,gen.k,scen.k,"absmid",trait.k] + R2selnAbsmid.pv
          #Constrain abs
          if(abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]>0.7) abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]=0.7
          if(abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]<0.4) abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]=0.4
          #rn evolution
          if(scen.k==5){
            abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]= abs.mean[yr.k,gen.k,scen.k,"rn",trait.k] + R2selnRN.pv
            #Constrain abs
            if(abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]>1) abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]=1
            if(abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]< -1) abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]=-1
          }
            
            #HB
            trait.k=2
            abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]= abs.mean[yr.k,gen.k,scen.k,"absmid",trait.k] + R2selnAbsmid.hb
            #Constrain abs
            if(abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]>0.7) abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]=0.7
            if(abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]<0.4) abs.mean[yr.k,gen.k+1,scen.k,"absmid",trait.k]=0.4
            #rn evolution
            if(scen.k==5){
              abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]= abs.mean[yr.k,gen.k,scen.k,"rn",trait.k] + R2selnRN.hb
              #Constrain abs
              if(abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]>1) abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]=1
              if(abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]< -1) abs.mean[yr.k,gen.k+1,scen.k,"rn",trait.k]=-1
          }
          
        } #end check generation
        
        #also put in next year's slot
        if(yr.k<length(years)-1){
          
        #PV
        trait.k=1  
        abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]= abs.mean[yr.k,gen.k,scen.k,"absmid",trait.k] + R2selnAbsmid.pv
        #Constrain abs
        if(abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]>0.7) abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]=0.7
        if(abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]<0.4) abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]=0.4
        
        if(scen.k==5) {abs.mean[yr.k+1,1,scen.k,"rn",trait.k]= abs.mean[yr.k,gen.k,scen.k,"rn",trait.k] + R2selnRN.pv
        #Constrain abs
        if(abs.mean[yr.k+1,1,scen.k,"rn",trait.k]>1) abs.mean[yr.k+1,1,scen.k,"rn",trait.k]=1
        if(abs.mean[yr.k+1,1,scen.k,"rn",trait.k]< -1) abs.mean[yr.k+1,1,scen.k,"rn",trait.k]=-1}
        
        #HB
        trait.k=2  
        abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]= abs.mean[yr.k,gen.k,scen.k,"absmid",trait.k] + R2selnAbsmid.hb
        #Constrain abs
        if(abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]>0.7) abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]=0.7
        if(abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]<0.4) abs.mean[yr.k+1,1,scen.k,"absmid",trait.k]=0.4
        
        if(scen.k==5) {abs.mean[yr.k+1,1,scen.k,"rn",trait.k]= abs.mean[yr.k,gen.k,scen.k,"rn",trait.k] + R2selnRN.hb
        #Constrain abs
        if(abs.mean[yr.k+1,1,scen.k,"rn",trait.k]>1) abs.mean[yr.k+1,1,scen.k,"rn",trait.k]=1
        if(abs.mean[yr.k+1,1,scen.k,"rn",trait.k]< -1) abs.mean[yr.k+1,1,scen.k,"rn",trait.k]=-1}
      
        } #end check year
        
        #Store other metrics
        abs.mean[yr.k,gen.k,scen.k,"abssample",1]= colMeans(abs.plast.pv)
        abs.mean[yr.k,gen.k,scen.k,"Babsmid",1]= BetaAbsmid.pv
        if(scen.k==5) abs.mean[yr.k,gen.k,scen.k,"Brn",1]= BetaRN.pv
        
        abs.mean[yr.k,gen.k,scen.k,"abssample",2]= colMeans(abs.plast.hb)
        abs.mean[yr.k,gen.k,scen.k,"Babsmid",2]= BetaAbsmid.hb
        if(scen.k==5) abs.mean[yr.k,gen.k,scen.k,"Brn",2]= BetaRN.hb
        
      } #end scen loop
      
  } #end generation
  print(yr.k)
} #end year 

#=====================================
#Save output
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/out/")
#saveRDS(abs.mean, "absmean.abs")
#saveRDS(lambda.mean, "lambdamean.abs")

#abs.mean <- readRDS("absmean.abs")
#lambda.mean <- readRDS("lambdamean.abs")

#Not much evolution

#lambda.mean[yr.k,gen.k,scen.k]
#abs.mean[yr.k,gen.k,scen.k,"abssample"]

#all scenarios
#pv
abs.scen= rbind(abs.mean[,,1,"abssample",1],abs.mean[,,2,"abssample",1],abs.mean[,,3,"abssample",1],abs.mean[,,4,"abssample",1],abs.mean[,,5,"abssample",1] )
abs.scen= as.data.frame(abs.scen)
abs.scen$scenario= c(rep("plast0evol0",33),rep("plast1evol0",33),rep("plast0evol1",33),rep("plast1evol1",33),rep("plast1evol1rnevol1",33) )
abs.scen$year= c(rep(years,5) )
colnames(abs.scen)[1:5]=c("gen1","gen2","gen3","gen4","gen5")
abs.scen$trait="pv"
abs.scen.pv= abs.scen
#hb
abs.scen= rbind(abs.mean[,,1,"abssample",2],abs.mean[,,2,"abssample",2],abs.mean[,,3,"abssample",2],abs.mean[,,4,"abssample",2],abs.mean[,,5,"abssample",2] )
abs.scen= as.data.frame(abs.scen)
abs.scen$scenario= c(rep("plast0evol0",33),rep("plast1evol0",33),rep("plast0evol1",33),rep("plast1evol1",33),rep("plast1evol1rnevol1",33) )
abs.scen$year= c(rep(years,5) )
colnames(abs.scen)[1:5]=c("gen1","gen2","gen3","gen4","gen5")
abs.scen$trait="hb"
abs.scen= rbind(abs.scen.pv,abs.scen)

#subset scenarios
abs.scen=abs.scen[which(abs.scen$scenario %in% c("plast0evol1","plast1evol0","plast1evol1") ),]
abs.scen$scen.name<- "evolution only"
abs.scen$scen.name[which(abs.scen$scenario=="plast1evol0" )]= "plasticity only"
abs.scen$scen.name[which(abs.scen$scenario=="plast1evol1" )]= "plasticity + evolution"
abs.scen$scen.name= factor(abs.scen$scen.name, levels=c("plasticity only","evolution only","plasticity + evolution") )

abs.l=melt(abs.scen, id.vars = c("year","scenario","scen.name","trait"))

fig.absevol= ggplot(abs.l, aes(x=year, y=value, color=variable, group=variable))+geom_line()+
  scale_color_viridis_d()+
  facet_grid(trait~scen.name)+
  theme_classic()+ylab("absorptivity (%)")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_AbsEvol.pdf", height = 4, width = 7)
fig.absevol
dev.off()

#just hb
fig.absevol= ggplot(abs.l[abs.l$trait=="hb",], aes(x=year, y=value, color=variable, group=variable))+geom_line()+
  scale_color_viridis_d()+
  facet_wrap(~scen.name)+
  theme_classic()+ylab("HB absorptivity (%)")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_AbsEvol_hb.pdf", height = 3, width = 7)
fig.absevol
dev.off()

#-----------
#pv
lambda.scen= rbind(lambda.mean[,,1,1],lambda.mean[,,2,1],lambda.mean[,,3,1],lambda.mean[,,4,1],lambda.mean[,,5,1] )
lambda.scen= as.data.frame(lambda.scen)
lambda.scen$scenario= c(rep("plast0evol0",33),rep("plast1evol0",33),rep("plast0evol1",33),rep("plast1evol1",33),rep("plast1evol1rnevol1",33) )
lambda.scen$year= c(rep(years,5) )
colnames(lambda.scen)[1:5]=c("gen1","gen2","gen3","gen4","gen5")
lambda.scen$trait="pv"
lambda.scen.pv= lambda.scen
#hb
lambda.scen= rbind(lambda.mean[,,1,2],lambda.mean[,,2,2],lambda.mean[,,3,2],lambda.mean[,,4,2],lambda.mean[,,5,2] )
lambda.scen= as.data.frame(lambda.scen)
lambda.scen$scenario= c(rep("plast0evol0",33),rep("plast1evol0",33),rep("plast0evol1",33),rep("plast1evol1",33),rep("plast1evol1rnevol1",33) )
lambda.scen$year= c(rep(years,5) )
colnames(lambda.scen)[1:5]=c("gen1","gen2","gen3","gen4","gen5")
lambda.scen$trait="hb"
lambda.scen= rbind(lambda.scen.pv, lambda.scen)

lambda.l=melt(lambda.scen, id.vars = c("year","scenario","trait"))

ggplot(lambda.l, aes(x=year, y=value, color=variable, group=variable))+geom_line()+
  scale_color_viridis_d()+
  facet_grid(trait~scenario)

