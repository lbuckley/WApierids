#R version 4.1.0 (2021-05-18)

library(ggplot2) #ggplot2_3.4.4
library(reshape2) #reshape2_1.4.4
library(reshape) #reshape_0.8.8
library(viridisLite) #viridisLite_0.4.0
library(patchwork) #patchwork_1.2.0
library(TrenchR) #TrenchR_1.1.1

#P. rapae larval TPC

#estimate larval temperature
#2003-2004
#2020-2021

#larval TPC
#growth rate
#Figure 2, mean short term mass-specific growth rate, from kingsolver 2000, 
#time period: July 25-30, Aug 12-19, Aug 21-26
#g/g/hr
temps= c(11,17,23,29,35,40,41)
mgr= c(.010, .018, .0435, .0494, .0726, .0388, .0165)
gr= as.data.frame(cbind(temps, mgr))

#selection gradient
#Figure 7
temps= c(11,17,23,29,35)
pm.sg= c(0.338,0.182, 0.190, -0.233, -0.217) #selection gradient for pupal mass
sg= as.data.frame(cbind(temps, pm.sg))
#add y values
sg$ys= gr$mgr[1:5]

#consumption rate
#Kingsolver 2000 PBZ
# 24h consumption rate for 4th instar, g/g/h
temps=c(10, 15, 20, 25, 30, 35, 40, 45)
cr= c(.012, .030, .0316, .0577, .0575, .077, .0428, .0085)
cr= as.data.frame(cbind(temps, cr))
#---------

#fit tpc

colnames(gr)=c("temp","rate")
d=gr

  #scale to max height
  d$rate_scale= d$rate*2.5/max(d$rate)
  
  fits <- nls(rate_scale~TPC_beta(T_b=temp, shift, breadth, tolerance, aran=0, skew),
              data = d,
              start= c(shift = 0, breadth = 0.1, tolerance=43, skew=0.6), 
              algorithm = "port",
              lower = c(shift = -9, breadth = 0.04, tolerance=20, skew=0.5),
              upper = c(shift = 9, breadth = 0.15, tolerance=60, skew=0.8)
  )
  
  tpc.beta= coef(fits)
  tpc.beta= c(tpc.beta, max(d$rate))
  
#plot with points
  #scale to max height
  d$rate_scale= d$rate*2.5/max(d$rate)
  
  plot(d$temp,d$rate_scale, ylim=c(0,3), xlim=c(0,50))
  points(1:50, TPC_beta(T_b=1:50, shift=tpc.beta[1], breadth=tpc.beta[2], aran=0, tolerance=tpc.beta[3], skew=tpc.beta[4]), type="l")
  
#--------
#plot P. rapae temperature distributions

locations= c("Corfu","Seattle")
loc.k=2

#years for data
if(loc.k==1) years=c(1989:2021) #1989:1993, 2017:2021
if(loc.k==2) years=c(1998:2021) #2001:2005, 2017:2021

#combine data
for(yr.k in 1:length(years)){
dat= read.csv(paste("./data/era5_micro_sun/",locations[loc.k],years[yr.k],".csv",sep="") )
dat$year= years[yr.k]

#vary periods
if(loc.k==1){
  if(years[yr.k]<2000)dat$period="initial"
  if(years[yr.k]>=2000 & years[yr.k]<2010)dat$period="middle"
  if(years[yr.k]>=2010)dat$period="recent"
}
if(loc.k==2){
  if(years[yr.k]<2006)dat$period="initial"
  if(years[yr.k]>=2006 & years[yr.k]<2015)dat$period="middle"
  if(years[yr.k]>=2015)dat$period="recent"
}

if(yr.k==1) dat.all=dat
if(yr.k>1) dat.all=rbind(dat, dat.all)
}

#subset to sunlight, keep 24h
dat.day= dat.all #subset(dat.all, dat.all$SOLR>0)

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
dat.day.plot= dat.day[-which(dat.day$period=="middle"),]

#--------------------------  
#plot density distributions
p1= ggplot(dat.day.plot, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Growth rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+
  theme(legend.position = c(0.1, 0.95))+
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 
#D0cm, TALOC, TAREF

#------------------
#Plot just for period of selection study: July 28-Aug 5, 9 days expand to 14 days, July 26-Aug 8 

dat.day1= dat.day[which(dat.day$DOY %in% c(207:220)),]

#relabel period
pers= c("initial","middle","recent")
yrs= c("1998-2005","2006-2013","2014-2021")

dat.day1$period= yrs[match(dat.day1$period, pers)]

#percent temps exceeding 34.5C
length(which(dat.day1$TALOC[dat.day1$period=="1998-2005"]>34.5))/length(dat.day1$TALOC)
length(which(dat.day1$TALOC[dat.day1$period=="2014-2021"]>34.5))/length(dat.day1$TALOC)
#percent temps exceeding 20C
length(which(dat.day1$TALOC[dat.day1$period=="1998-2005"]>20))/length(dat.day1$TALOC)
length(which(dat.day1$TALOC[dat.day1$period=="2014-2021"]>20))/length(dat.day1$TALOC)

#drop middle period
dat.day1.plot= dat.day1[-which(dat.day1$period=="2006-2013"),]

#plot density distributions
p1= ggplot(dat.day1.plot, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  ylab("Density")+
  xlab("Temperature (°C)" )+
  theme_classic(base_size = 20) +
  theme(legend.position = "none")+ #theme(legend.position = c(0.2, 0.9))+
  xlim(5,50)+
  scale_y_continuous(expand=c(0,0))+
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 

p1.ref= ggplot(dat.day1.plot, aes(x=TAREF))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  ylab("Density")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20) +
  theme(legend.position = c(0.15, 0.9))+
  xlim(5,45)+
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 

#plot together
p1.pl.ref=p1 + geom_density(alpha=0.4, linetype="dashed", aes(x=TAREF, color=period))+
  theme(legend.position = c(0.7, 0.9))+
  ylim(0,0.095)

#===============================
#P. rapae larvae
#plot TPC
p.dat= as.data.frame(cbind(1:50, TPC_beta(1:50, shift=tpc.beta[1], breadth=tpc.beta[2], aran=0, tolerance=tpc.beta[3], skew=tpc.beta[4])/(2.5/tpc.beta[5])  ))
colnames(p.dat)= c("Tb","performance")
p.dat$performance[is.nan(p.dat$performance)]=0

p2= p1 + geom_line(data=p.dat, aes(x = Tb, y = performance) )

# add selection arrows
p3= p2 + geom_segment(data=sg, aes(x = temps, y = ys, xend = temps, yend = ys+pm.sg/20),
                      arrow = arrow(length = unit(0.3, "cm")), lwd=1)+
  ylim(0,0.095)

#add points 
p3= p3 + geom_point(data=d, aes(x = temp, y = rate))+
  theme(legend.position = c(0.7, 0.98))

#analogous plot for reference height
p3.ref= p1.ref + geom_line(data=p.dat, aes(x = Tb, y = performance) )+ 
  geom_segment(data=sg, aes(x = temps, y = ys, xend = temps, yend = ys+pm.sg/20),
  arrow = arrow(length = unit(0.3, "cm")), lwd=1)+
  geom_point(data=d, aes(x = temp, y = rate))

#-----------
#just study period
dat.day=dat.day1

#estimate performance
dat.day$perf= TPC_beta(dat.day$TALOC, shift=tpc.beta[1], breadth=tpc.beta[2], aran=0, tolerance=tpc.beta[3], skew=tpc.beta[4])/(2.5/tpc.beta[5])

#plot density distributions of performance
p1= ggplot(dat.day, aes(x=perf))+
  geom_density(alpha=0.5, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  xlab("Performance")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.75, 0.95))+
  scale_color_manual(values=c("#7AD151FF","#440154FF")  )+
  scale_fill_manual(values=c("#7AD151FF","#440154FF")  ) 

#Count of NAs above CTmax of TPC
#tab= table( is.nan(dat.day$perf), dat.day$period)
#prop.CTmax= c( #tab["TRUE","initial"]/(tab["FALSE","initial"]+tab["TRUE","initial"]),
#               tab["TRUE","middle"]/(tab["FALSE","middle"]+tab["TRUE","middle"]),
#               tab["TRUE","recent"]/(tab["FALSE","recent"]+tab["TRUE","recent"]))

#performance means
dat.day1= na.omit(dat.day)
aggregate(dat.day1$perf, list(dat.day1$period), FUN=mean)

#density distribution by performance
dat.day1$Tround= round(dat.day1$TALOC)
dat.day2= aggregate(dat.day1$perf, list(dat.day1$Tround, dat.day1$period, dat.day1$seas), FUN=sum)
names(dat.day2)= c("temp","period","seas","sumperf")
dat.day2$per_seas= paste(dat.day2$period, dat.day2$seas, sep="_")

#normalize to count of data
counts= aggregate(dat.day1$perf, list(dat.day1$period, dat.day1$seas), FUN=function(x)length(x))
names(counts)= c("period","seas","count")
counts$per_seas= paste(counts$period, counts$seas, sep="_")
#add counts to performance data
dat.day2$counts= counts$count[match(dat.day2$per_seas, counts$per_seas)]
#normalize
dat.day2$sumperf.norm= dat.day2$sumperf/dat.day2$counts

p4= ggplot(dat.day2, aes(x=temp, y=sumperf.norm, color=period))+ #geom_line()+   
  geom_smooth(se=FALSE)+
 # facet_wrap(~seas)+
  xlab("Temperature at reference height (°C)")+
  ylab("Sum of feeding rate (g/g/h)" )+
  #ylim(0,0.0015)+xlim(0,40)+
  theme_classic(base_size = 20)+theme(legend.position = c(0.4, 0.2))

pdf("./figures/Fig_Prapae_FeedingRateByTemp.pdf", height = 6, width = 6)
p4
dev.off()

#-----------------
#P. rapae
#Kingsolver JG (2000) Feeding, growth, and the thermal environment of cabbage white caterpillars, Pieris rapae L. Physiological and Biochemical Zoology, 73(5):621–628.

#http://labs.bio.unc.edu/Buckley/WGdocs/Kingsolver2001.pdf
#selection on pupal mass, development time, survival

#TPCs for growth and development rates will be used to estimate size and development time as well as the # microclimates experienced by pupae and adults as a result of developmental timing 

#estimate growth and development with different TPCs to estimate selection (vary beta parameters?)
#estimate for each season and year to estimate selection through time

params= expand.grid(temp = seq(0,70,1), shift= seq(tpc.beta[1]-10,tpc.beta[1]+10,1),
                    breadth=c(tpc.beta[2]-0.03,tpc.beta[2],tpc.beta[2]+0.03) )
  
beta.mat= function(pmat,tolerance, skew) TPC_beta(pmat[1], pmat[2], pmat[3], aran=0, tolerance, skew)

params$perf= apply(params, MARGIN=1, FUN=beta.mat, tolerance=tpc.beta[3], skew=tpc.beta[4] )
params$perf[is.nan(params$perf)]= 0

#---

#estimate growth and development
#RUN
param.grid= expand.grid(shift= seq(tpc.beta[1]-15,tpc.beta[1]+10,1),
                        breadth=c(tpc.beta[2]-0.03,tpc.beta[2],tpc.beta[2]+0.03) )

perf.mat= matrix(NA, nrow= nrow(dat.day), ncol= nrow(param.grid) )
for(topt.k in 1:nrow(param.grid)){
  perf.mat[,topt.k]= sapply(dat.day$TALOC, FUN=TPC_beta, 
         shift=param.grid[topt.k,1], breadth=param.grid[topt.k,2], aran=0, 
         tolerance=tpc.beta[3], skew=tpc.beta[4])
  #set NaN to zero
  perf.mat[which(is.nan(perf.mat[,topt.k])),topt.k]= 0
  #account for scaling
  perf.mat[,topt.k]= perf.mat[,topt.k]/(2.5/tpc.beta[5])
}

#estimate topt for curve
param.grid$Topt=NA
for(topt.k in 1:nrow(param.grid)){
  ps= TPC_beta(1:50,shift=param.grid[topt.k,1], breadth=param.grid[topt.k,2], aran=0, 
               tolerance=tpc.beta[3], skew=tpc.beta[4])
  param.grid$Topt[topt.k]= c(1:50)[which.max(ps)]
}

#combine dates
perfs= cbind(dat.day[,c("year","period","seas")], perf.mat)

#aggregate
#not seasons, just study period
perfs1= aggregate(perfs[,4:ncol(perfs)], list(perfs$year,perfs$period), FUN=mean)
names(perfs1)=c("year","period",1:nrow(param.grid))

#plot fitness curves through time
perfs.l= melt(perfs1[,1:(ncol(perfs)-1)], id.vars = c("year","period"))
names(perfs.l)[3:4]=c("temperature","performance")

perfs.l$shift= param.grid$shift[perfs.l$temperature]
perfs.l$breadth= param.grid$breadth[perfs.l$temperature]
perfs.l$topt= param.grid$Topt[perfs.l$temperature]

#perfs.l= melt(perfs1[,1:ncol(perfs1)], id.vars = c("year","period"))
#names(perfs.l)[3:4]=c("temperature","performance")
perfs.l$temperature= as.numeric(as.character(perfs.l$temperature))

#group by breadth
perfs.l$yrbr= paste(perfs.l$year, perfs.l$breadth,sep="_")
#make breadth a factor
perfs.l$breadth= factor(perfs.l$breadth)

#combine shift, breadth
fig.fitnesscurves=ggplot(perfs.l[which(perfs.l$breadth==0.15),], aes(x=topt, y=performance, color=year, group=yrbr) )+geom_line()+
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+
  xlab("Thermal optima (C)")+ylab("Mean feeding rate (g/g/h)") +
  theme(legend.position = c(0.2, 0.8))+
  xlim(5,45)

fig.fitnesscurves.all=ggplot(perfs.l, aes(x=topt, y=performance, color=year, group=yrbr) )+geom_line(aes(lty=breadth))+
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+
  xlab("thermal optima (C)")+ylab("feeding rate (g/g/h)") +
  xlim(5,45)

#over time
fig.pyr=ggplot(perfs.l[which(perfs.l$breadth==0.15),], aes(x=year, y=performance, color=topt, group=topt) )+geom_line()+
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+
  xlab("Thermal optima (C)")+ylab("Mean feeding rate (g/g/h)") +
  theme(legend.position = c(0.7, 0.8))+
  xlim(5,45)

#find thermal optima through years  
breadths= sort(unique(param.grid$breadth))

inds= which(param.grid$breadth==breadths[1])+3
perfs1$shift_opt_b1= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

inds= which(param.grid$breadth==breadths[2])+3
perfs1$shift_opt_b2= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

inds= which(param.grid$breadth==breadths[3])+3
perfs1$shift_opt_b3= param.grid$Topt[apply(perfs1[,inds], MARGIN=1, FUN=which.max )]

#melt
perfs.b= perfs1[,c("year","period","shift_opt_b1","shift_opt_b2","shift_opt_b3")]
perfs.b.l= melt(perfs.b, id.vars = c("year","period"))
names(perfs.b.l)[3:4]=c("breadth","opt_shift")
perfs.b.l$breadths= NA
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b1"]<- breadths[1]
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b2"]<- breadths[2]
perfs.b.l$breadths[perfs.b.l$breadth=="shift_opt_b3"]<- breadths[3]

perfs.b.l$seas_br= paste(perfs.b.l$seas, perfs.b.l$breadths,sep="_")
perfs.b.l$breadth= factor(perfs.b.l$breadths)

fig.shift_opt= ggplot(perfs.b.l[which(perfs.b.l$breadth==0.15),], aes(x=year, y=opt_shift, group= seas_br))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  xlab("Year")+ylab("Optimal Topt (C)")+
  guides(lty="none")+ 
  scale_color_manual(values=c("#39568CFF", "#DCE319FF"))

fig.shift_opt.all= ggplot(perfs.b.l, aes(x=year, y=opt_shift, group= seas_br, lty=breadth ))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  xlab("Year")+ylab("Optimal Topt (C)")+ 
  scale_color_manual(values=c("#39568CFF", "#DCE319FF"))

#Figure 5
pdf("./figures/Fig5_PrapaeStudy.pdf", height = 12, width = 8)
p3 / fig.fitnesscurves / fig.shift_opt +
  plot_layout(heights = c(2, 1.5, 1.25))+
  plot_annotation(tag_levels = 'A')
dev.off()

pdf("./figures/Fig_PrapaeStudy_supp.pdf", height = 10, width = 8)
fig.fitnesscurves.all / fig.shift_opt.all +
  plot_layout(heights = c(1.5, 1.25))+
  plot_annotation(tag_levels = 'A')
dev.off()

temps.wa<- p1.pl.ref

#----------------------------
##P. rapae 1999 selection data
# rgr at each temperature
# survival to pupation, time to pupation, pupal mass
#survive to eclosion, time to eclosion, adult mass
#number of eggs laid
#Relate to selection estimates: determine Topt and plot against fitness components, only pupal mass was significant?

#compare to empirical data
pr= read.csv("./data/KingsolverPrapae/PrapaeUW.Seln2.1999.Combineddata.OPUS2021.csv")

#plot TPCs
pr1= pr[,c("Mom","UniID", "Mi", "RGR11", "RGR17", "RGR23", "RGR29", "RGR35")]
pr1= melt(pr1, id.vars=c("Mom","UniID", "Mi"), variable.name="temp", value.name="rgr")           
pr1$temperature= gsub('RGR', '', pr1$temp)

ggplot(pr1, aes(x=temperature, y=rgr, color=UniID, group=UniID))+geom_line()+
  facet_wrap(~Mom)

#fit TPC
#load FitGaussian

fitG =
  function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }

ids= unique(pr1$UniID)
#by mom
#ids= unique(pr1$Mom)

tpc.p= matrix(NA, nrow=length(ids), ncol=3)

for(id.k in 1:length(ids)){

#gr= pr1[pr1$Mom==ids[id.k],c("temperature","rgr")]
gr= pr1[pr1$UniID==ids[id.k],c("temperature","rgr")]
colnames(gr)=c("temp","rate")
gr$temp= as.numeric(gr$temp)
gr= na.omit(gr)

#gr.fit= try( fitG(x=gr$temp, y=gr$rate, mu=31, sig=8, scale= 0.06) )
#tpc.p[id.k,]=gr.fit$par

tryCatch({ gr.fit<- fit.tpcs(gr) 
tpc.p[id.k,]<-coef(gr.fit) },
error=function(e){})

}
#fix so Topt can be >35

#plot TPCs
pr$Topt= tpc.p[match(pr$UniID,ids), 2]

#plot selection functions
ggplot(pr, aes(x=Topt, y=Time.to.Pupation, color=UniID, group=UniID))+geom_point()
ggplot(pr, aes(x=Topt, y=Pupa.wt, color=UniID, group=UniID))+geom_point()
ggplot(pr, aes(x=Topt, y=Butt..Wt, color=UniID, group=UniID))+geom_point()
ggplot(pr, aes(x=Topt, y=Fecundity, color=UniID, group=UniID))+geom_point()

#============================
#PLOT PARAMETER COMBINATIONS

tpc1<- TPC_beta(T_b=1:50, shift=0, breadth=0.15, aran=0, tolerance=30, skew=0.70)

tpc2<- TPC_beta(T_b=1:50, shift=5, breadth=0.15, aran=0, tolerance=30, skew=0.70)

tpc3<- TPC_beta(T_b=1:50, shift=0, breadth=0.25, aran=0, tolerance=30, skew=0.70)

tpc4<- TPC_beta(T_b=1:50, shift=0, breadth=0.15, aran=0, tolerance=40, skew=0.70)

names(tpc1)[2]<- "value"
names(tpc2)[2]<- "value"
names(tpc3)[2]<- "value"
names(tpc4)[2]<- "value"

tpc1<- as.data.frame(cbind(1:50, tpc1))
tpc1$parameters<- "Tmin= 0, breadth=0.15, b=30"
tpc2<- as.data.frame(cbind(1:50, tpc2))
tpc2$parameters<- "Tmin= 5, breadth=0.15, b=30"
tpc3<- as.data.frame(cbind(1:50, tpc3))
tpc3$parameters<- "Tmin= 0, breadth=0.25, b=30"
tpc4<- as.data.frame(cbind(1:50, tpc4))
tpc4$parameters<- "Tmin= 0, breadth=0.15, b=40"

names(tpc1)[2]<- "value"
names(tpc2)[2]<- "value"
names(tpc3)[2]<- "value"
names(tpc4)[2]<- "value"

tpcs<- rbind(tpc1, tpc2, tpc3, tpc4)

names(tpcs)[1]<- "temp"

tpc.param<- ggplot(tpcs, aes(x=temp, y=value, color=parameters))+
  geom_line(lwd=1.5)+
  theme_classic(base_size=20)+ scale_color_viridis_d()+
  xlab("Body temperature (C)")+ylab("Performance")+
  theme(legend.position = "right") #+guides(fill=guide_legend(nrow=2,byrow=TRUE))

#setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("./figures/FigS_TPC.pdf", height = 6, width = 10)
tpc.param
dev.off()

