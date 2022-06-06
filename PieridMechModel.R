## TO DO
# Add size into fitness estimates 
#convert from wing traits to absorptivity

library(truncnorm)
library(dplyr)
library(TrenchR)
library(ggplot2)
library(MCMCglmm)

#LOAD PARAMETERS
#Demographic parameters
#Kingsolver 1983 Ecology 64
OviRate=0.73; # Ovipositing rates: 0.73 eggs/min (Stanton 1980) 
MaxEggs=700; # Max egg production: 700 eggs (Tabashnik 1980)
PropFlight= 0.5; # Females spend 50# of available activity time for oviposition-related
# Watt 1979 Oecologia
SurvDaily=0.6; # Daily loss rate for Colias p. eriphyle at Crested Butte, female values
#Hayes 1981, Data for Colias alexandra at RMBL
SurvMat=0.014; #1.4# survival to maturity

#read/make species data
solar.abs= 0.65 # Solar absorptivity, proportion
SpecDat= as.data.frame(solar.abs)
SpecDat$d=0.36 # d- Thoractic diameter, cm
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
if(loc.k==2) years=c(2001:2005, 2017:2021)

setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro/')
#combine data
for(yr.k in 1:length(years)){
  dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
  dat$year= years[yr.k]
  if(years[yr.k]<2000)dat$period="initial"
  if(years[yr.k]>=2000 & years[yr.k]<2010)dat$period="middle"
  if(years[yr.k]>=2010)dat$period="recent"
  
  if(yr.k==1) dat.all=dat
  if(yr.k>1) dat.all=rbind(dat, dat.all)
}

#subset to sunlight
dat.day= subset(dat.all, dat.all$TIME>260)
dat.day= subset(dat.all, dat.all$TIME<1080)

#divide by periods, 91 to 274
#April to May
dat.day$seas="AprMay"
#June to July
dat.day$seas[dat.day$DOY>151]="JunJul"
#August to September
dat.day$seas[dat.day$DOY>212]="AugSep"
#order
dat.day$seas= factor(dat.day$seas, levels=c("AprMay","JunJul","AugSep") )

#combine doy and time
dat.day$d.hr= dat.day$DOY + dat.day$TIME/60

#subset columns
dat.sub= dat.day[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","ZEN","SOLR","VLOC","D0cm","year","period","seas")]
#TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
#TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
#VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)

#melt temp data
datm= melt(dat.sub[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","D0cm","year","period")], id=c("dates","DOY","TIME","d.hr","year","period") )

#plot
ggplot(data=datm, aes(x=d.hr, y = value, color=variable))+ geom_line(alpha=0.2)+
  theme_bw()+scale_color_viridis_d()

#partition solar radiation [or extract from micro?], returns diffuse fraction
df=partition_solar_radiation("Erbs", kt=0.8)

#butterfly temperature
dat.sub$Tb= Tb_butterfly( T_a = dat.sub$TALOC, Tg = dat.sub$D0cm, Tg_sh = dat.sub$D0cm, u = dat.sub$VLOC, 
                          H_sdir = dat.sub$SOLR*(1-df), H_sdif = dat.sub$SOLR*(df), z = dat.sub$ZEN, D = 0.36, 
                          delta = 1.46, alpha = 0.6, r_g = 0.3, wing_angle=42)
#fix high
dat.sub$Tb[which(dat.sub$Tb>80)]<-NA

#plot Tb distributions
ggplot(data=dat.sub, aes(x=d.hr, y = Tb))+ geom_line(alpha=0.4)+
  theme_bw()+scale_color_viridis_d()
#plot density distributions
p1= ggplot(dat.sub, aes(x=Tb))+
  geom_density(alpha=0.5, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Body Temperature (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.6, 0.8))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccTb.pdf", height = 7, width = 10)
p1
dev.off()

#---------
### SET UP DATA STRUCTURES
#Make array to store data
#Lambda= survival*fecundity
#Matrix of lambdas
#dims: stats, year, Lambda 
Lambda<-array(NA, dim=c(length(years),length(seq(0.4,0.7,0.05)),5,5)) #Last dimension is Lambda, FAT,Egg Viability, growth
dimnames(Lambda)[[1]]<-years
dimnames(Lambda)[[3]]<-c("gen1","gen2","gen3","gen4","gen5")

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
 
  # ESTIMATE DEVELOPMENTAL TIMING
  #https://www.jstor.org/stable/4913
  
  DevZeros= c(10) #egg to pupal
  GddReqs= c(105) 
  
  dat.yr$dd= dat.yr$TALOC -DevZeros
  dat.yr$dd[dat.yr$dd<0]= 0 #set negative to zero
  dat.yr$dd= dat.yr$dd/24 #divide by 24 to swirch from degree days to hours
  dat.yr$dd.cs= cumsum(dat.yr$dd)
  
  #estiamte pupal timing
  #compare to empirical data
  setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/')
  pr= read.csv("PrapaeUW.Seln2.1999.Combineddata.OPUS2021.csv")
  
  #mean time to pupation and eclosion
  colMeans(pr[,c("Time.to.Pupation","Time.to.Eclosion")], na.rm=TRUE)
  #about 7 day pupation
  
  #-----------------------------------
  #Calculate development timing
  
    for(gen.k in 1:5){
      
      #Assume 7 days from eclosion to eggs laid
      #Hatching ~5days 
      Jlarv= ifelse(gen.k>1, Jadult+12, dat.yr$DOY[which.max(dat.yr$dd>0)] )  
      if(Jlarv>max(dat.yr$DOY)) Jlarv=max(dat.yr$DOY)
      
      ##TO PUPATION
      inds= which(dat.yr$DOY==Jlarv)[1]:nrow(dat.yr)
      Jpup= dat.yr$DOY[inds[which.max(cumsum(dat.yr$dd[inds])> GddReqs)]] 
       if(Jpup>max(dat.yr$DOY) | length(Jpup)==0) Jpup=max(dat.yr$DOY) 
      
      #PUPATION 
      # ~7 day pupation basedd on CUH field data
      Jadult= Jpup + 7
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
  facet_wrap(~variable)
plot.pocc.temps= ggplot(pup.temps.l[which(pup.temps.l$variable %in% c("Tlarv","Tpup","Tad") ),], aes(x=yr, y=value, col=gen, group=gen))+geom_line()+
  facet_wrap(~variable)

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
  
  for(gen.k in 1:5 ){ #loop generation
    
      for(abs.k in 1:length(abs1) ){ #loop absorptivity
        
        #ADD PLASTICITY
        #https://doi.org/10.1093/icb/38.3.545
        #Basal dorsal HW melanism
        #use females, check units of basal dorsal HW melanism
        
        #hours light
        light.h= c(16,10)
        abs.m= c(.66,.67); abs.m= c(.717,.735)        
        mod1= lm(abs.m~light.h)
        abs.plast= function(light.hrs, abs) abs -0.003*light.hrs
        
        light.hrs=daylength(lat=46.8151, doy=pup.temps["Jpup", yr.k, gen.k])
        
        #absorptivity
        abs= abs.plast(light.hrs, abs1[abs.k])
        
        #get flight dates
        Jfl=pup.temps["Jadult",yr.k, gen.k]   
        
        #calculate Tb based on absorptivity
        dat.yr$Tb.a= Tb_butterfly( T_a = dat.yr$TALOC, Tg = dat.yr$D0cm, Tg_sh = dat.yr$D0cm, 
                                  u = dat.yr$VLOC, 
                                  H_sdir = dat.yr$SOLR*(1-df), H_sdif = dat.yr$SOLR*(df), 
                                  z = dat.yr$ZEN, D = 0.36, 
                                  delta = 1.46, alpha = abs, r_g = 0.3)
        #fix high
        dat.yr$Tb.a[which(dat.yr$Tb.a>80)]<-NA
        
        #Flight probability
        dat.yr$fl.p= sapply(dat.yr$Tb.a, FUN=fl.ph)
        
        #Egg viability
        dat.yr$egg.v= sapply(dat.yr$Tb.a, FUN=egg.viab)
      
        #larval growth
        #Use P. rapae curve
        tpc.gaus= c(0.06340297, 30.92347862,  8.40085316)
        dat.yr$larv.growth= sapply(dat.yr$TALOC, FUN=gaussian_1987, rmax=tpc.gaus[1], topt=tpc.gaus[2], a=tpc.gaus[3])
        
       # daily values
        dat.day<- dat.yr %>%
          group_by(DOY) %>%
          summarise(FAT = sum(fl.p), EggViab= geo_mean(egg.v), Tb.mean= mean(Tb.a), Gt.l = sum(larv.growth) )
        
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
        
        if(!is.nan(Eggs)){
          MaxDay=5
          Lambda1=0
          for(day in 1:MaxDay){
            Eggs1= min(Eggs, MaxEggs-Eggs*(day-1))  ###LIMIT MAX NUMBER EGGS
            if(Eggs1<0) Eggs1=0
            Lambda1= Lambda1+ SurvMat * SurvDaily^day *Eggs1;                        
          }#end loop days
          
          Lambda[yr.k, abs.k, gen.k, ]= c(Lambda1, mean(FAT.ind), mean(ev.ind), mean(T.ind, na.rm=T), Gr.larv )
          
        } #Check Eggs
        
      } #end loop absorptivity
     
  } #end loop generation
  
} #end loop years

#SAVE OBJECT
#filename= paste("lambda1_",projs[proj.k],".rds",sep="")
#saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#melt lambda array
Lambda.l= melt(Lambda) 
colnames(Lambda.l)=c("year","abs","gen","metric","value")
Lambda.l$abs= abs1[Lambda.l$abs]
Lambda.l$metric= c("lambda","fat","ev","tind","gr.l")[Lambda.l$metric]

#plot lambdas
Lambda.l$gyr= paste(Lambda.l$gen, Lambda.l$year, sep="")

fig.OccFit= ggplot(Lambda.l, aes(x=abs, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccFitCurv.pdf", height = 7, width = 10)
fig.OccFit
dev.off()

#-------------
#test Tb function
dat.sub1= dat.sub[dat.sub$DOY==100 & dat.sub$year==2021,]
Ts= matrix(NA, nrow=nrow(dat.sub1), ncol=length(abs1))

for(abs.k in 1:length(abs1) ){
  Ts[,abs.k]= Tb_butterfly( T_a = dat.sub1$TALOC, Tg = dat.sub1$D0cm, Tg_sh = dat.sub1$D0cm, 
                              u = dat.sub1$VLOC, 
                              H_sdir = dat.sub1$SOLR*(1-df), H_sdif = dat.sub1$SOLR*(df), 
                              z = dat.sub1$ZEN, D = 0.36, 
                              delta = 1.46, alpha = abs1[abs.k], r_g = 0.3, wing_angle=42, shade = FALSE) 
}

Ts= as.data.frame(Ts)
colnames(Ts)= abs1

Ts$hr= dat.sub1$TIME/60+12
Ts$hr[Ts$hr>24]= Ts$hr[Ts$hr>24]-24
Ts$TALOC= dat.sub1$TALOC

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

abs.opt= array(NA, dim=c(length(years), 5))  

lambda.fit=function(x, a=a) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients
lambda.max=function(x, a=a) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] 

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    Lambda.yr.gen= Lambda[yr.k, , gen.k, ]
    
    #Extract temperatures
    Tp= pup.temps["Tpup",yr.k, gen.k]
    
    #--------------------------
    #Fitness models
    
    if(!all(is.na(Lambda.yr.gen[,1]))){ #check has data
      
      fit= lm(Lambda.yr.gen[,1]~a+I(a^2))
      #Estimate fitness functions across cells
      fit.mod= fit$coefficients
      
      #find maxima lambda
      abs.opt[yr.k,gen.k]= a.fit$a[which.max(predict.lm(fit, a.fit))]
      
    } #end check data
  } #end gen loop
} #end year loop

#save optimal Alphas
#saveRDS(abs.opt, paste("abs.opt_",projs[proj.k],".rds", sep=""))
#abs.opt <- readRDS( paste("abs.opt_",projs[proj.k],".rds", sep="") )

#plot optima alphas
abs.opt.l=as.data.frame(abs.opt)
abs.opt.l$year=years
abs.opt.l=melt(abs.opt.l,id.vars = c("year"))

ggplot(abs.opt.l, aes(x=year, y=value, color=variable, group=variable))+geom_line()+
  scale_color_viridis_d()

#***************************************
#compute initial AbsMean 
#UPDATE
int_elev = 0.4226; slope_elev = 0.06517
Tmid = 20; slope_plast = -0.0083  #if Tmid=22.5, -0.006667;

h2= 0.4
#abs.mean=0.55
abs.sd= 0.062
rn.sd= 0.0083

## NEED TO CALC ABS.OPT
#initialize with optimum value across for ten years, across generations
abs.init2 <- mean(abs.opt[1:10, ], na.rm=TRUE)
#Use optimal
abs.init<- abs.init2

#-----------------------
#Save values
abs.mean= array(NA, dim=c(length(years), 5, 5,5))  #dims: yr.k, gen.k, scen.k:no plast, plast, only plast, metrics: abssample, absmid, rn, Babsmid, Brn)
abs.mean[1,,,2]= abs.init
abs.mean[1,,,3]= slope_plast
dimnames(abs.mean)[[4]]= c("abssample", "absmid", "rn", "Babsmid", "Brn") 

lambda.mean= array(NA, dim=c(length(years),5, 5)) #dims: yr.k, gen.k, scen.k:no plast, plast, only plast)

#BetaRN= rep(NA, nrow(pts.sel))
#-------------------------------
scen.mat= rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,0),c(1,1,1) )
colnames(scen.mat)= c("plast","evol","evolRN"  )

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    BetaAbsmid=NA
    
    Lambda.yr.gen= Lambda[yr.k, , gen.k, ]
    
    #determine those completing generations
    comp.gen= which(pup.temps["Jadult",yr.k,gen.k]<243)
    nocomp.gen= which(pup.temps["Jadult",yr.k,gen.k]==243)
    #set those not completing generations to NA
    if(length(nocomp.gen)>0) Lambda.yr.gen[nocomp.gen,]=NA
    
      #Extract temperatures
      Tp= pup.temps["Tpup",yr.k, gen.k]
      
      #Fitness models
      fit= lm(Lambda.yr.gen[,1]~a+I(a^2))
      #Estimate fitness functions across cells
      fit.mod= fit$coefficients
      
      #find maxima lambda
      abs.max= a.fit$a[which.max(predict.lm(fit, a.fit))]
      
      #-------------------------
      # LOOP PLASTICITY SCENARIOS
      for(scen.k in 1:5){ #plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1
        
        if(scen.mat[scen.k,1]==1) rn.mean1= slope_plast
        if(scen.mat[scen.k,1]==0) rn.mean1= 0
        if(scen.k==5 & gen.k==1) rn.mean1= abs.mean[yr.k,gen.k,scen.k,"rn"]
        if(scen.k==5 & gen.k>1) rn.mean1= abs.mean[yr.k,gen.k-1,scen.k,"rn"]
        
        if(gen.k==1) abs.mean1= abs.mean[yr.k,gen.k,scen.k,"absmid"]
        if(gen.k>1) abs.mean1= abs.mean[yr.k,gen.k-1,scen.k,"absmid"]
        
        #change NA values to negative values 
        abs.na.inds= abs.mean1[which( is.na(abs.mean1))]
        rn.na.inds= rn.mean1[which( is.na(rn.mean1))]
        
        abs.mean1[which( is.na(abs.mean1))]= -10 
        rn.mean1[which( is.na(rn.mean1))]= -1000
        
        #Choose random sample of abs and rn values from current distribution (truncated normal) 
        abs.sample= sapply(abs.mean1, function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
        rn.sample= sapply(rn.mean1, function(x) rtnorm(N.ind, mean = x, sd = rn.sd, lower=-1, upper=1) )
        if(scen.mat[scen.k,1]==0) rn.sample[]=0
        
        #Add plasticity across sites and sample
        abs.plast <- abs.sample + rn.sample*(Tp-Tmid)
        #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
        
        ##calculate fitness 
        #use fitness function to predict Lambda for each individual
        #extract coefficients and calculate across abs samples
        fit.sample=  sapply(abs.plast[,1], function(x) fit.mod[1]+x*fit.mod[2]+x^2*fit.mod[3] )
       
        #Fit.pred <- eval.fd(Abs.sample,Fitmod.year.gen) ### for spline
        
        #standardize to relative fitness and centered on trait mean
        fit.mean= mean(fit.sample)
        lambda.mean[yr.k,gen.k,scen.k]=fit.mean
        rel.fit= fit.sample/fit.mean
        
        absmid.dif= t( apply(abs.sample,1,'-',abs.mean1) )
        rn.dif= t( apply(rn.sample,1,'-',rn.mean1) )
        
        R2selnAbsmid<- 0 #No response to selection if no evolution
        R2selnRN<- 0 
        #------------
        if(scen.k<5 & scen.mat[scen.k,2]==1){
          ##selection analysis
          rel.fit= as.vector(rel.fit)
          absmid.dif= as.vector(absmid.dif)
          sel.fit= lm(rel.fit~absmid.dif +I(absmid.dif^2))$coefficients
          sel.mod= lm(rel.fit~absmid.dif +I(absmid.dif^2))
          
          #Response to selection
          BetaAbsmid <-sel.fit[2]
          R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
        } #end scen.k<5
        #------------
        if(scen.k==5){    
          ##selection analysis
          rel.fit= as.vector(rel.fit)
          absmid.dif= as.vector(absmid.dif)
          sel.fit= lm(rel.fit~absmid.dif +I(absmid.dif^2))$coefficients
          sel.mod= lm(rel.fit~absmid.dif +I(absmid.dif^2))
          
          #Response to selection
          BetaAbsmid <-sel.fit[2]
          R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
          
          BetaRN <- sel.fit[3] 
          R2selnRN <- h2*(rn.sd^2)*BetaRN
        } #end scen.k==5
        #-------------
        
        #Response to selection
        if(gen.k<5) {
          abs.mean[yr.k,gen.k+1,scen.k,"absmid"]= abs.mean[yr.k,gen.k,scen.k,"absmid"] + R2selnAbsmid
          #Constrain abs
          if(abs.mean[yr.k,gen.k+1,scen.k,"absmid"]>0.7) abs.mean[yr.k,gen.k+1,scen.k,"absmid"]=0.7
          if(abs.mean[yr.k,gen.k+1,scen.k,"absmid"]<0.4) abs.mean[yr.k,gen.k+1,scen.k,"absmid"]=0.4
          
          #rn evolution
          if(scen.k==5){
            abs.mean[yr.k,gen.k+1,scen.k,"rn"]= abs.mean[yr.k,gen.k,scen.k,"rn"] + R2selnRN
            #Constrain abs
            if(abs.mean[yr.k,gen.k+1,scen.k,"rn"]>1) abs.mean[yr.k,gen.k+1,scen.k,"rn"]=1
            if(abs.mean[yr.k,gen.k+1,scen.k,"rn"]< -1) abs.mean[yr.k,gen.k+1,scen.k,"rn"]=-1
          }
        } #end evolutionary scenarios
        
        #also put in next year's slot
        if(yr.k<length(years)-1){
        abs.mean[yr.k+1,1,scen.k,"absmid"]= abs.mean[yr.k,gen.k,scen.k,"absmid"] + R2selnAbsmid
        #Constrain abs
        if(abs.mean[yr.k+1,1,scen.k,"absmid"]>0.7) abs.mean[yr.k+1,1,scen.k,"absmid"]=0.7
        if(abs.mean[yr.k+1,1,scen.k,"absmid"]<0.4) abs.mean[yr.k+1,1,scen.k,"absmid"]=0.4
        
        if(scen.k==5) {abs.mean[yr.k+1,1,scen.k,"rn"]= abs.mean[yr.k,gen.k,scen.k,"rn"] + R2selnRN
        #Constrain abs
        if(abs.mean[yr.k+1,1,scen.k,"rn"]>1) abs.mean[yr.k+1,1,scen.k,"rn"]=1
        if(abs.mean[yr.k+1,1,scen.k,"rn"]< -1) abs.mean[yr.k+1,1,scen.k,"rn"]=-1}
        }
        
        #Store other metrics
        abs.mean[yr.k,gen.k,scen.k,"abssample"]= colMeans(abs.plast)
        abs.mean[yr.k,gen.k,scen.k,"Babsmid"]= BetaAbsmid
        if(scen.k==5) abs.mean[yr.k,gen.k,scen.k,"Brn"]= BetaRN
        
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
abs.scen= rbind(abs.mean[,,1,"abssample"],abs.mean[,,2,"abssample"],abs.mean[,,3,"abssample"],abs.mean[,,4,"abssample"],abs.mean[,,5,"abssample"] )
abs.scen= as.data.frame(abs.scen)
abs.scen$scenario= c(rep("plast0evol0",33),rep("plast1evol0",33),rep("plast0evol1",33),rep("plast1evol1",33),rep("plast1evol1rnevol1",33) )
abs.scen$year= c(rep(years,5) )
colnames(abs.scen)[1:5]=c("gen1","gen2","gen3","gen4","gen5")

abs.l=melt(abs.scen, id.vars = c("year","scenario"))

ggplot(abs.l, aes(x=year, y=value, color=variable, group=variable))+geom_line()+
  scale_color_viridis_d()+
  facet_wrap(~scenario)

#-----------

lambda.scen= rbind(lambda.mean[,,1],lambda.mean[,,2],lambda.mean[,,3],lambda.mean[,,4],lambda.mean[,,5] )
lambda.scen= as.data.frame(lambda.scen)
lambda.scen$scenario= c(rep("plast0evol0",33),rep("plast1evol0",33),rep("plast0evol1",33),rep("plast1evol1",33),rep("plast1evol1rnevol1",33) )
lambda.scen$year= c(rep(years,5) )
colnames(lambda.scen)[1:5]=c("gen1","gen2","gen3","gen4","gen5")

lambda.l=melt(lambda.scen, id.vars = c("year","scenario"))

ggplot(lambda.l, aes(x=year, y=value, color=variable, group=variable))+geom_line()+
  scale_color_viridis_d()+
  facet_wrap(~scenario)

