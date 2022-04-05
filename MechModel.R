library(ggplot2)
library(reshape2)
library(reshape)

# Adapt Colias niche model

#global historical climatology network?

#Center for urban horticulture
#lat 47.657628, Longitude: -122.290255

#Corfu site in central Washington in the Columbia National Wildlife Refuge
#12 miles west of Othello
#46.850663264 -119.535331192
  
#------------------


#-----------------------------------------------
#P. rapae larval TPC

#estimate larval temperature
#2003-2004
#2020-2021

#larval TPC
#growth rate
#Figure 2, mean short term mass-specific growth rate, from kingsolver 2000, 
#mg/g/hr
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

# plot
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

#estimate performance shift

#model selection

#--------
#plot temperature distributions

locations= c("Corfu","Seattle")
loc.k=1

#years for data
if(loc.k==1) years=c(1989:1993, 2017:2021)
if(loc.k==2) years=c(2001:2005, 2017:2021)

setwd('/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/data/era5_micro/')
#combine data
for(yr.k in 1:10){
dat= read.csv(paste(locations[loc.k],years[yr.k],".csv",sep="") )
dat$year= years[yr.k]
if(yr.k<6)dat$period="initial"
if(yr.k>5)dat$period="recent"

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
  
#plot time series
p1= ggplot(dat.day, aes(x=d.hr, y=TALOC, color=year))+
  geom_line(alpha=0.2)+
  xlab("Temperature at reference height (°C)")+
  ylab("Consumption or growth rate (g/h/h)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.75, 0.9))

#plot density distributions
p1= ggplot(dat.day, aes(x=TALOC))+
  geom_density(alpha=0.5, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  xlab("Temperature at reference height (°C)")+
  ylab("Consumption or growth rate (g/h/h)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.75, 0.9))
#D0cm, TALOC, TAREF
#Why is TAREF consistently lower than TALOC? Over water?

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

#estimate performance
dat.day$perf= beta_2012(dat.day$TALOC, tpc.beta[1], tpc.beta[2], tpc.beta[3], tpc.beta[4], tpc.beta[5])

#plot density distributions
p1= ggplot(dat.day, aes(x=perf))+
  geom_density(alpha=0.5, aes(fill=period, color=period))+
  facet_wrap(~seas)+
  xlab("Performance")+
  ylab("Density" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.75, 0.9))

#Count of NAs above CTmax of TPC
table( is.nan(dat.day$perf), dat.day$period)

#============================
#P. occidentalis adult selection

#subset columns
dat.sub= dat.day[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","ZEN","SOLR","VLOC","D0cm")]
#TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
#TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
#VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)

#melt temp data
datm= melt(dat.sub[,c("dates","DOY","TIME","d.hr","TAREF","TALOC","D0cm")], id=c("dates","DOY","TIME","d.hr") )

#plot
ggplot(data=datm, aes(x=d.hr, y = value, color=variable))+ geom_line(alpha=0.2)+
  theme_bw()+ 
  #xlab("ordinal date") +ylab("abundance")+ 
  #labs(color = "seasonal GDDs")+ theme(strip.text = element_text(face = "italic")) 
  scale_color_viridis()

#butterfly temperature

#partition solar radiation [or extract from micro?]










