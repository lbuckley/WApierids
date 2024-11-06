#Simplify model to correspond to survey

#Kingsolver JG (1995) Viability selection on seasonally polyphenic traits: wing melanin pattern in western white butterflies. Evolution, 49(5):932–941

#new figure 6. seasonal variation and selection
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/data/")
sv= read.csv("DF_Kingsolver1995b_Fig2.csv")
betas= read.csv("Survival_Kingsolver1995b_Table3-4.csv")

plot.sv= ggplot(sv, aes(x=doy, y=Discriminant.score, shape=Sex))+
  geom_point(size=4) + theme_classic(base_size = 18) +
  geom_line()+
  theme(legend.position = c(0.2,0.35))+
  ylab("Discriminant Function Score")+
  xlab("Day of Year")+
  scale_shape_manual(values=c(1,2))

#Add selection arrows

#make offsets for different traits
offs=c(-3,0,3)
betas$pos= betas$doy + offs[match(betas$Trait, c("pv","hb","fwl") )] 
betas$Trait= factor(betas$Trait, levels=c("pv", "hb","fwl") )

#ys
doys= c(130, 140, 183, 205)
y.doys= c(0.663, 0.663,-1.6, -1.98)
betas$ys= y.doys[match(betas$doy, doys )] 

#plot
plot.sv2= plot.sv + geom_segment(data=betas, aes(x = pos, y = ys, xend = pos, 
                                                 yend = ys+1.5*selection.coefficient, 
                                                 color=Trait, shape="f"),
                                 arrow = arrow(length = unit(0.3, "cm")), lwd=1.2)+
  scale_color_manual(values=alpha(c("#66c2a5", "#fc8d62","#8da0cb"),1) )
#scale_color_manual(values=c("#440154FF", "#238A8DFF","#B8DE29FF") )

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig7_seas.pdf",height = 5, width = 4)
plot.sv2
dev.off()

#study periods
#two week periods
#June 21-28, 1989 (N =289 marked butterflies); 169-182
#June 27-July 7, 1990 (N = 247); 177-190
#May 15-24, 1991 (N = 403); 133-146
#May 5-15, 1992  (N = 147) 123-136

#------------------------------
#plot temperature distributions
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

#subset columns
dat.sub= dat.day[,c("dates","DOY","TIME","TAREF","TALOC","ZEN","SOLR","VLOC","D0cm","year","period")]
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

#subset columns
dat.sub.sh= dat.day.sh[,c("dates","DOY","TIME","TAREF","TALOC","ZEN","SOLR","VLOC","D0cm","year","period")]
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

#----------
#butterfly temperature
Temat= cbind(dat.sub$TALOC, dat.sub$TALOC_sh,dat.sub$D0cm, dat.sub$D0cm_sh, 
             dat.sub$VLOC, dat.sub$SOLR*(1-df), dat.sub$SOLR*(df), dat.sub$ZEN)

dat.sub$Tb= apply(Temat, MARGIN=1, FUN=Tb_butterfly.mat,
                  D = 0.26, delta = 1.46, HB = 0.55, PV = 0.55, 
                  r_g = 0.3, wing_angle=42, shade = TRUE) 

#fix high
dat.sub$Tb[which(dat.sub$Tb>80)]<-NA

##add study periods
#dat.sub$DOYs=NA
#dat.sub$DOYs[dat.sub$DOY %in% c(169:182)]="June 18-July 2"
#dat.sub$DOYs[dat.sub$DOY %in% c(177:190)]="June 26-July 9"
#dat.sub$DOYs[dat.sub$DOY %in% c(133:146)]="May 13-26"
#dat.sub$DOYs[dat.sub$DOY %in% c(123:136)]="May 3-16"

#combine spring and summer periods, 4 weeks
dat.sub$seas=NA
dat.sub$seas[dat.sub$DOY %in% c(169:192)]="June-July"
dat.sub$seas[dat.sub$DOY %in% c(123:146)]="May"
dat.sub$seas= factor(dat.sub$seas, levels=c("May","June-July") )

#restrict to study periods
dat.sub= dat.sub[which(!is.na(dat.sub$seas)),]

#plot density distributions
p1= ggplot(dat.sub, aes(x=Tb))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Body Temperature (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 16)+theme(legend.position = c(0.45, 0.7))+
  xlim(0,90)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccTb.study.pdf", height = 7, width = 10)
p1
dev.off()

#air temp
p2= ggplot(dat.sub, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Environmental Temperature (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.4, 0.85))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccTaloc.study.pdf", height = 7, width = 10)
p2
dev.off()

#========================
#adults of given size

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
  
  print(years[yr.k])
  
  dat.yr= subset(dat.sub, dat.sub$year==years[yr.k])
  
  Temat= cbind(dat.yr$TALOC, dat.yr$TALOC_sh,dat.yr$D0cm, dat.yr$D0cm_sh, 
               dat.yr$VLOC, dat.yr$SOLR*(1-df), dat.yr$SOLR*(df), dat.yr$ZEN)
  
  for(gen.k in 1:2 ){ #loop generation
    
    for(abs.k.pv in 1:length(abs1) ){ #loop PV absorptivity
      for(abs.k.hb in 1:length(abs1) ){ #loop HB absorptivity
        
        abs.pv= abs1[abs.k.pv]
        abs.hb= abs1[abs.k.hb]
        
        #calculate Tb based on absorptivity
        dat.yr$Tb.a= apply(Temat, MARGIN=1, FUN=Tb_butterfly.mat,
                           D = 0.26, delta = 1.46, HB = abs.hb, PV = abs.pv, 
                           r_g = 0.3, wing_angle=42, shade = TRUE) 
        
        #Flight probability
        dat.yr$fl.p= sapply(dat.yr$Tb.a, FUN=fl.ph)
        
        #Egg viability
        dat.yr$egg.v= sapply(dat.yr$Tb.a, FUN=egg.viab)
        
        # daily values
        dat.day<- dat.yr %>%
          group_by(DOY) %>%
          summarise(FAT = sum(fl.p), EggViab= geo_mean(egg.v), Tb.mean= mean(Tb.a),
                    .groups="keep")
        
        #flight day
        Jfl= c(128,174)[gen.k]
        
        ##CALCULATE EGG VIABILITY OVER 10 DAY PERIOD (GEOMETRIC MEAN ACROSS HOURS)
        #sample flight day from truncated normal distribution
        Nind=200 #1000 #changed from 100
        f.low= Jfl+5
        f.up= Jfl-5
        
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
        
        Eggs.ind= 60*PropFlight*OviRate*FAT.ind * ev.ind #account for Egg viability
        Eggs.ind_noViab= 60*PropFlight*OviRate*FAT.ind
        
        #Means across individuals
        Eggs= mean(Eggs.ind)
        Eggs_noViab= mean(Eggs.ind_noViab)
        
        #Set max eggs as a function of ppupal size
        #https://doi.org/10.1071/ZO9820223
        #Fig 5 in Jones et al. life time fecundity, y= 4.84w -273, w=weight(mg)
        #bound pupal weight 110-210 mg
        Gr.larv=160 #assume mean value
        MaxEggs= 4.84*Gr.larv-273
        
        #Survival to Maturity
        #egg surivival 0.18, Davies and Gilbert 1985, https://www.jstor.org/stable/pdf/4217728.pdf
        #assume 20 days larval duration
        #P. virginiensis larval survival, https://doi.org/10.1890/05-0647
        SurvMat= 0.18*larv.surv(20)
        
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
#filename= paste("lambda1_",projs[proj.k],".rds",sep="")
#saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#melt lambda array
Lambda.l= melt(Lambda) 
colnames(Lambda.l)=c("year","abs.pv","abs.hb","gen","metric","value")
Lambda.l$abs.pv= abs1[Lambda.l$abs.pv]
Lambda.l$abs.hb= abs1[Lambda.l$abs.hb]
Lambda.l$gen= c("May","June-July")[match(as.character(Lambda.l$gen), c("gen1","gen2"))]
Lambda.l$metric= c("fitness","FAT (h)","egg viab (%)","Tadult (C)","larval growth")[Lambda.l$metric]
Lambda.l$metric= factor(Lambda.l$metric, levels=c("fitness","FAT (h)","egg viab (%)","Tadult (C)","larval growth"))

#plot lambdas
Lambda.l$gyr= paste(Lambda.l$gen, Lambda.l$year, sep="")

#first two gens
Lambda.l= Lambda.l[which(!is.na(Lambda.l$gen)),]
#order generations
Lambda.l$gen= factor(Lambda.l$gen, levels=c("May","June-July") )

fig.OccFit.pv= ggplot(Lambda.l[which(Lambda.l$abs.hb==0.55 & Lambda.l$metric %in% c("fitness","FAT (h)","egg viab (%)")),], 
                      aes(x=abs.pv, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()+
  theme_classic()+xlab("PV ventral absorptivity (%)")+ylab("")

fig.OccFit.hb= ggplot(Lambda.l[which(Lambda.l$abs.pv==0.55 & Lambda.l$metric %in% c("fitness","FAT (h)","egg viab (%)")),], 
                      aes(x=abs.hb, y=value, color=year, group=gyr))+geom_line()+
  facet_grid(metric~gen, scale="free_y") +scale_color_viridis_c()+
  theme_classic(base_size = 16)+xlab("HB dorsal absorptivity (%)")+ylab("")+
  theme(legend.position="bottom",legend.key.width=unit(0.1,"npc"))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccFitCurv.pv.study.pdf", height = 7, width = 7)
fig.OccFit.pv
dev.off()

pdf("Fig_PoccFitCurv.hb.study.pdf", height = 10, width = 6)
p1 + fig.OccFit.hb +
  plot_layout(heights = c(0.8, 2))+ 
  plot_annotation(tag_levels = 'A')
dev.off()

#==============================================
#Estimate selection gradients

ngens=2
sel.abs= array(NA, dim=c(length(years), 2,2)) #last dimension is pv, hb  

lambda.B=function(x, a=abs1) if(sum(is.na(x))==0) lm(x~a)$coefficients

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, 1]
    
    #to long format
    Lambda.fit= melt(Lambda.yr.gen) 
    colnames(Lambda.fit)=c("abs.pv","abs.hb","lambda")
    Lambda.fit$abs.pv= abs1[Lambda.fit$abs.pv]
    Lambda.fit$abs.hb= abs1[Lambda.fit$abs.hb]
    
    sel.abs[yr.k, gen.k,]= lm(Lambda.fit$lambda~Lambda.fit$abs.pv+Lambda.fit$abs.hb)$coefficients[2:3]
      
  } #end gen loop
} #end year loop

#plot
sel.abs.l= melt(sel.abs)
colnames(sel.abs.l)= c("years","gen","abs.metric","B")
sel.abs.l$years= years[sel.abs.l$years]
sel.abs.l$gen= c("May","June-July")[sel.abs.l$gen]
sel.abs.l$abs.metric= c("PV ventral melanism","HB dorsal melanism")[sel.abs.l$abs.metric]

sel.abs.l$gen= factor(sel.abs.l$gen,levels= c("May","June-July"))

#plot selection gradient through time
fig.B= ggplot(sel.abs.l, aes(x=years, y=B))+geom_line()+
  facet_grid(abs.metric~gen)+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  ylab("selection gradient") +xlab("year")+
  geom_hline(yintercept=0, color = "gray")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PoccB.study.pdf", height = 7, width = 7)
fig.B
dev.off()







