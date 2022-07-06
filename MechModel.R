library(ggplot2)
library(reshape2)
library(reshape)
library(viridisLite)
library(patchwork)

# Adapt Colias niche model
# https://github.com/lbuckley/ColiasBiogeog/blob/master/ColiasMain.R

#global historical climatology network?

#Center for urban horticulture
#lat 47.657628, Longitude: -122.290255

#Corfu site in central Washington in the Columbia National Wildlife Refuge
#12 miles west of Othello
#46.850663264 -119.535331192
  
#-----------------------------------------------
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
#rTPC, https://github.com/padpadpadpad/rTPC
#https://padpadpadpad.github.io/rTPC/articles/fit_many_models.html

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#fit absolute growth rate
#https://doi.org/10.1111/j.0014-3820.2004.tb01732.x
#Fig 2, mg/hr for 4th instaar
#temps= c(11,17,23,29,35,40)
#growth= c(0.192,0.448,0.993,1.111,1.461,0.700)
#gr= as.data.frame(cbind(temps, growth))

colnames(gr)=c("temp","rate")
d=gr

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
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:weibull)

# get parameters using tidy
params <- d_stack %>%
  mutate(., est = map(fit, tidy)) %>%
  select(-fit) %>%
  unnest(est)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
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

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

#extract coefficients
tpc.beta= coef(mod)
#beta_2012(temp, a, b, c, d, e)

#extract weibull model
#mod= d_fits$weibull[[1]]
#coef(mod)

#--------
#plot P. rapae temperature distributions

locations= c("Corfu","Seattle")
loc.k=2

#years for data
if(loc.k==1) years=c(1989:2021) #1989:1993, 2017:2021
if(loc.k==2) years=c(1998:2021) #2001:2005, 2017:2021

#SUN
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_sun/')
#combine data
for(yr.k in 1:length(years)){
dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
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

#subset to sunlight
dat.day= subset(dat.all, dat.all$SOLR>0)

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

dat.day.sun= dat.day

#----------------
#SHADE
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro_shade/')
#combine data
for(yr.k in 1:length(years)){
  dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
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

#subset to sunlight
dat.day= subset(dat.all, dat.all$SOLR>0)

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

dat.day.sh= dat.day

#--------------------------  
#plot density distributions
p1= ggplot(dat.day.sun, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  ylab("Growth rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.65, 0.8))
#D0cm, TALOC, TAREF

#use sun
dat.day= dat.day.sun

#------------------
#Plot just for period of selection study: July 28-Aug 5, 9 days expand to 14 days, July 26-Aug 8 

dat.day1= dat.day[which(dat.day$DOY %in% c(207:220)),]

#relabel period
pers= c("initial","middle","recent")
yrs= c("1998-2005","2006-2013","2014-2021")

dat.day1$period= yrs[match(dat.day1$period, pers)]

#plot density distributions
p1= ggplot(dat.day1, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  ylab("Growth rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.22, 0.9))

## get 1999 (study) data
#plot just during study
ggplot(dat.day1[dat.day1$year==1999,], aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+
  ylab("Growth rate (g/g/h)")+
  xlab("Temperature at plant height (°C)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.2, 0.8))

#===============================
#P. rapae larvae
#plot TPC
p.dat= as.data.frame(cbind(1:50, beta_2012(1:50, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])))
colnames(p.dat)= c("Tb","performance")
p.dat$performance[is.nan(p.dat$performance)]=0

p2= p1 + geom_line(data=p.dat, aes(x = Tb, y = performance) )

# add selection arrows
#i + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
#                 arrow = arrow(length = unit(0.5, "cm")))
p3= p2 + geom_segment(data=sg, aes(x = temps, y = ys, xend = temps, yend = ys+pm.sg/20),
                      arrow = arrow(length = unit(0.3, "cm")), lwd=1)

#-----------
#just study period
dat.day=dat.day1

#estimate performance
dat.day$perf= beta_2012(dat.day$TALOC, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])

#plot density distributions of performance
p1= ggplot(dat.day, aes(x=perf))+
  geom_density(alpha=0.5, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  xlab("Performance")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.75, 0.9))

#Count of NAs above CTmax of TPC
tab= table( is.nan(dat.day$perf), dat.day$period)
prop.CTmax= c( #tab["TRUE","initial"]/(tab["FALSE","initial"]+tab["TRUE","initial"]),
               tab["TRUE","middle"]/(tab["FALSE","middle"]+tab["TRUE","middle"]),
               tab["TRUE","recent"]/(tab["FALSE","recent"]+tab["TRUE","recent"]))

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
  geom_smooth()+
 # facet_wrap(~seas)+
  xlab("Temperature at reference height (°C)")+
  ylab("Sum of growth rate (g/g/h)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.7, 0.8))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_Prapae.pdf", height = 10, width = 12)
p3 / p4
dev.off()

#-----------------
#P. rapae
#Kingsolver JG (2000) Feeding, growth, and the thermal environment of cabbage white caterpillars, Pieris rapae L. Physiological and Biochemical Zoology, 73(5):621–628.

#http://labs.bio.unc.edu/Buckley/WGdocs/Kingsolver2001.pdf
#selection on pupal mass, development time, survival

#TPCs for growth and development rates will be used to estimate size and development time as well as the # microclimates experienced by pupae and adults as a result of developmental timing 

#estimate growth and development with different TPCs to estimate selection (vary beta parameters?)
#estimate for each season and year to estimate selection through time

#vary TPC
#Gaussian, vary Topt
mod= d_fits$gaussian[[1]]
tpc.gaus= coef(mod)

#generate parameter combinations
temps=1:70

params= expand.grid(temp = seq(0,70,2), topt = seq(10, 30, 1) )
  
gauss.mat= function(pmat,rmax,a) gaussian_1987(pmat[1], rmax, pmat[2], a)

params$perf= apply(params, MARGIN=1, FUN=gauss.mat, rmax=tpc.gaus[1], a=tpc.gaus[3])

ggplot(params, aes(x=temp, y=perf, color=topt, group=topt))+geom_line()

#estimate growth and development
#RUN
topts= seq(10, 30, 1)
perf.mat= matrix(NA, nrow= nrow(dat.day), ncol= length(topts) )
for(topt.k in 1:length(topts)){
  perf.mat[,topt.k]= sapply(dat.day$TALOC, FUN=gaussian_1987, rmax=tpc.gaus[1], topt=topts[topt.k], a=tpc.gaus[3])
}

#combine dates
perfs= cbind(dat.day[,c("year","period","seas")], perf.mat)
#aggregate
#perfs1= aggregate(perfs[,4:24], list(perfs$year,perfs$period,perfs$seas), FUN=mean)
#names(perfs1)=c("year","period","seas", seq(10, 30, 1))

#not seasons, just study period
perfs1= aggregate(perfs[,4:24], list(perfs$year,perfs$period), FUN=mean)
names(perfs1)=c("year","period",seq(10, 30, 1))

#make column for slopes
perfs1$B= NA

for(row.k in 1: nrow(perfs1)){
  #put into format for regression
  perfs2= as.data.frame(cbind(10:30, t(perfs1[row.k,4:24]) ))
  colnames(perfs2)=c("topt","perf")
  
  plot(perfs2$topt, perfs2$perf)
  #mod1= lm(perf~topt , data=perfs2)
  mod1= lm(perf~topt +I(topt^2) , data=perfs2)
  perfs1$B[row.k]= coefficients(mod1)[2]
  perfs1$curve[row.k]= coefficients(mod1)[3]
} 

#plot fitness curves through time
#perfs.l= melt(perfs1[,1:24], id.vars = c("year","period","seas"))
#names(perfs.l)[4:5]=c("temperature","performance")

perfs.l= melt(perfs1[,1:24], id.vars = c("year","period"))
names(perfs.l)[3:4]=c("temperature","performance")
perfs.l$temperature= as.numeric(as.character(perfs.l$temperature))

fig.fitnesscurves=ggplot(perfs.l, aes(x=temperature, y=performance, color=year, group=year))+geom_line()+
  #facet_wrap(~seas) +
  scale_color_viridis_c()+
  theme_classic(base_size = 20)+
  xlab("thermal optima (C)")+ylab("growth rate (g/g/h)")+theme(legend.position = c(0.6, 0.3))

#plot selection gradients through time
fig.selb=ggplot(perfs1, aes(x=year, y=B, color=seas))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(se=FALSE)
fig.selcurve= ggplot(perfs1, aes(x=year, y=curve, color=seas))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(se=FALSE)

#drop season
fig.selb=ggplot(perfs1, aes(x=year, y=B))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(se=FALSE)
fig.selcurve= ggplot(perfs1, aes(x=year, y=curve))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(se=FALSE)

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PrapaeSelection.pdf", height = 12, width = 10)
fig.fitnesscurves / fig.selb / fig.selcurve
dev.off()

#account for temperatures exceeding tpcs

#find thermal optima through years 
perfs1$Topt= apply(perfs1[,4:26], MARGIN=1, FUN=which.max )

ggplot(perfs1, aes(x=year, y=Topt, color=seas))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)

#drop season
perfs1$Topt= apply(perfs1[,3:25], MARGIN=1, FUN=which.max )

fig.topt= ggplot(perfs1, aes(x=year, y=Topt))+geom_line()+
  theme_classic(base_size = 20)+geom_smooth(method="lm",se=FALSE)+
  ylab("thermal optima (C)")

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_PrapaeStudy.pdf", height = 8, width = 14)
p3 + fig.fitnesscurves+fig.topt +
  plot_layout(widths = c(2, 1.5, 1.25))
dev.off()

#----------------------------
##P. rapae 1999 selection data
# rgr at each temperature
# survival to pupation, time to pupation, pupal mass
#survive to eclosion, time to eclosion, adult mass
#number of eggs laid
#Relate to selection estimates: determine Topt and plot against fitness components, only pupal mass was significant?

#compare to empirical data
setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/')
pr= read.csv("PrapaeUW.Seln2.1999.Combineddata.OPUS2021.csv")

#plot TPCs
pr1= pr[,c("Mom","UniID", "Mi", "RGR11", "RGR17", "RGR23", "RGR29", "RGR35")]
pr1= melt(pr1, id.vars=c("Mom","UniID", "Mi"), variable.name="temp", value.name="rgr")           
pr1$temperature= gsub('RGR', '', pr1$variable)

ggplot(pr1, aes(x=temperature, y=value, color=UniID, group=UniID))+geom_line()+
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

#gr= pr1[pr1$Mom==ids[id.k],c("temperature","value")]
gr= pr1[pr1$UniID==ids[id.k],c("temperature","value")]
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

#P. occidentalis adult selection

#SEE PieridMechModel.R

#estimate selection on survival and fecundity from mechanistic model
#Fig 6: pv, hb, fwl

#Kingsolver JG (1995) Viability selection on seasonally polyphenic traits: wing melanin pattern in western white butterflies. Evolution, 49(5):932–941.
#direction of selection on one wing trait important for thermoregulation, melanin on the base of the dorsal hindwings (trait hb), fluctuated seasonally; there was evidence of directional selection for increased hb in the spring studies and for decreased hb in the summer studies

#Butterfly model
#estimate larval TPCs: development, size

#tpupal -> plasticity
#_> Tadult

#->Adult performance: survival and fecundity from mechansitic model
#how to relate to field data?

##Kingsolver 1995a
#natural collection
#5 wing melanism traits and FWL -> survival
##Kingsolver 1995b
#wing melanism traits (?) and FWL -> survival

#-----
#Development
#P. occidentalis RMBL, only 1 temp, 2 temps P. napi
#https://www.jstor.org/stable/pdf/4215136.pdf

#P. napi devlopment temp
#https://www.jstor.org/stable/pdf/3544871.pdf

#P. rapae vancouver
#** https://doi.org/10.2307/4537, https://www.jstor.org/stable/4913
#other populations: https://www.journals.uchicago.edu/doi/full/10.1086/317758
#egg to pupae 10C, 105 dd
#https://doi.org/10.1111/j.1420-9101.2007.01318.x

#size and egg production
#https://www.publish.csiro.au/ZO/ZO9820223

#Phylogeny
#https://academic.oup.com/biolinnean/article/88/3/413/2691608

#Plasticity: photoperiod plasticity
# Kingsolver and Wiernasz. 1991. Seasonal polyphenism. https://doi.org/10.1086/285195

#Mech model parameters
#P. napi flights as a function of Tb, Kingsolver 1985. Thermal ecology. https://www.jstor.org/stable/pdf/4217668.pdf











