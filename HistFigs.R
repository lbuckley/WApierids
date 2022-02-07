library(ggplot2)
library(patchwork)
library(viridis)

#Fig 5. P.rapae: TPC, operative temperatures, arrows indicating selection on pupal mass. Mostly selection at low temps
#Consumption rate
#4th and 5th instar growth rates

#operative temperatures
#calendar date 223 to 237, Seattle, WA, 28 Jul to 5 Aug 1999
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/")
to= read.csv("UWCUH.MetData.Seln2.Aug1999.csv", na.strings="-6999")
#average models
to$Tmodel= rowMeans(to[,c(6:20)], na.rm=TRUE)

#selection on growth
#data extracted from graph
#http://labs.bio.unc.edu/Buckley/WGdocs/Kingsolver2001.pdf

#Figure 2, mean short term mass-specific growth rate, from kingsolver 2000, 
#mg/g/hr
temps= c(11,17,23,29,35,40,41)
mgr= c(.010, .018, .0435, .0494, .0726, .0388, .0165)
gr= as.data.frame(cbind(temps, mgr))

#Figure 7
temps= c(11,17,23,29,35)
pm.sg= c(0.338,0.182, 0.190, -0.233, -0.217) #selection gradient for pupal mass
sg= as.data.frame(cbind(temps, pm.sg))

#---------
#Kingsolver 2000 PBZ
# 24h consumption rate for 4th instar, g/g/h
temps=c(10, 15, 20, 25, 30, 35, 40, 45)
cr= c(.012, .030, .0316, .0577, .0575, .077, .0428, .0085)
cr= as.data.frame(cbind(temps, cr))

#24hr growth rate
temps=c(10, 15, 20, 25, 30, 35, 40)
rgr4= c(.00976, .0145, .0324, .041, .0478, .0513, .0277)
rgr5= c(.0070, .0158, .0311, .050, .0536, .0617, .0042)
gr2a= as.data.frame(cbind(temps, rgr4))
gr2a$instar="4th"
names(gr2a)[2]="rgr"
gr2b= as.data.frame(cbind(temps, rgr5))
gr2b$instar="5th"
names(gr2b)[2]="rgr"
gr2= rbind(gr2a, gr2b)

#---------
#Kingsolver and Gomulkiewicz. 2003. Environmental variation and selection. https://doi.org/10.1093/icb/43.3.470
#Consumption rate 24h

#Figure 1
#6hr growth rate for 4th and 5th instar
#4th and 5th instar
temps= c(11,17,23,29,35,40)
rgr4= c(0.004, 0.0072, 0.0195, 0.0224, 0.0325, 0.0174)
gr3= as.data.frame(cbind(temps, rgr4))

temps= c(11,17,23,29,35,41)
rgr5= c(0.00535, 0.011, 0.0176, 0.022, 0.0263, 0.0015)
gr4= as.data.frame(cbind(temps, rgr5))
#---------------------
#density plot
p1= ggplot(to, aes(x=Tmodel))+
  geom_density(fill="gray", color="gray")+
  xlim(0,45)+
  xlab("Temperature (°C)")+
  ylab("Consumption or growth rate (g/h/h)" )+
  theme_classic(base_size = 20)+theme(legend.position = c(0.6, 0.7))

# add selection arrows
#i + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
#                 arrow = arrow(length = unit(0.5, "cm")))
p2= p1 + geom_segment(data=sg, aes(x = temps, y = .09, xend = temps, yend = .09+pm.sg/20),
                      arrow = arrow(length = unit(0.5, "cm")))

#add growth rate
p3= p2 + geom_line(data=gr, aes(x=temps, y= mgr), color="#482677FF", lwd=1.2)

#add consumption
plot.tpc= p3 + geom_line(data=cr, aes(x=temps, y= cr), color="#29AF7FFF", lwd=1.2)+
  annotate("text",label="consumption rate ", 
           x = 34, y = 0.11, size = 6, colour = "#29AF7FFF")+
  annotate("text",label="growth rate ", 
           x = 34, y = 0.1, size = 6, colour = "#482677FF")

#check other TPCS
p5= plot.tpc + geom_line(data=gr2, aes(x=temps, y= rgr, col=instar), lwd=1.2)
p6= p5 + geom_line(data=gr3, aes(x=temps, y= rgr4), color="orange", lwd=1.2)
p7= p6 + geom_line(data=gr4, aes(x=temps, y= rgr5), color="purple", lwd=1.2)

#----------------------
#plot out
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/figures/")
pdf("TPCplot.pdf",height = 6, width = 8)
plot.tpc
dev.off()

#====================
#new figure 6. seasonal variation and selection
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/data/")
sv= read.csv("DF_Kingsolver1995b_Fig2.csv")
betas= read.csv("Survival_Kingsolver1995b_Table3-4.csv")

plot.sv= ggplot(sv, aes(x=doy, y=Discriminant.score, shape=Sex))+
  geom_point(size=3.5) + theme_classic(base_size = 18) +
  geom_line()+
  theme(legend.position = c(0.2,0.35))+
  ylab("Discriminant Function Score")+
  xlab("Day of Year")

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
                      arrow = arrow(length = unit(0.3, "cm")))+
  scale_color_manual(values=alpha(c("#66c2a5", "#fc8d62","#8da0cb"),1) )
  #scale_color_manual(values=c("#440154FF", "#238A8DFF","#B8DE29FF") )

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/figures/")
pdf("Fig7_seas.pdf",height = 5, width = 4)
plot.sv2
dev.off()

#===================
#Figure 7. Field selection figure

#Kingsolver. Am Nat 1996. Experimental manipulation. https://doi.org/10.1086/285852
#2nd time period throughout
#sharpie

#Kingsolver. 1995. Fitness consequences.  https://doi.org/10.1111/j.1558-5646.1995.tb02329.x
#photoperiod

#Table 1: is discriminant factor difference detween LD and SD?
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/historical/data/")
wt= read.csv("WingTraits_Kingsolver1995_Table1.csv")

#subset to summer
wt= wt[wt$season=="summer",]

#subset to significant traits
wt= wt[wt$trait %in% c("fwl","hb","pv"),]

#change trait names
wt$Trait<- "FW length"
wt$Trait[wt$trait=="hb"]<- "BDHW"
wt$Trait[wt$trait=="pv"]<- "VHW"

#make year factor
wt$year<- as.factor(wt$year)
#order traits
wt$Trait<- factor(wt$Trait, levels=c("VHW","BDHW", "FW length"), ordered=TRUE)

plot.wt= ggplot(wt, aes(x=Trait, y=value, col=year, shape=sex))+
  geom_point(size=3.5) + theme_classic(base_size = 18) +
  theme(legend.position ="none")+
  ylab("DF loadings for short-day vs long-day")+
  xlab("Wing Trait")+
  scale_color_manual(values=c("#440154FF", "#238A8DFF") )+
  scale_x_discrete(expand=c(0.05,0.2))
  
#Or Fig 2? with discriminant distribution

#Fig 3-5 behavior
beh= read.csv("Behavior_Kingsolver1995_Fig3-5.csv")
beh$group= paste(beh$year, beh$timeperiod, beh$sex, beh$behavior,sep="_")

#restrict to flight 
#also basking?
beh= beh[beh$behavior %in% c("flight"),] #"basking",

#make year factor
beh$year<- as.factor(beh$year)
#recode photoperiod
beh$Photoperiod= "long-day"
beh$Photoperiod[beh$photoperiod=="SD"]= "short-day"

plot.beh= ggplot(beh, aes(x=Photoperiod, y=fraction, col=year, shape=sex, group=group, lty=timeperiod))+
  geom_point(size=3.5) + geom_line(lwd=1.5)+
  #facet_grid(.~behavior)+
  theme_classic(base_size = 18)+
  ylab("Flight proportion") +xlab("Photoperiod")+
  scale_color_manual(values=c("#440154FF", "#238A8DFF") )+ 
  scale_x_discrete(expand=c(0.05,0.05))+
  scale_linetype_discrete(name="time of day")

#Fig 6-8 survival
#photoperiod
s.p= read.csv("Survival_Kingsolver1995_Fig6-8.csv")
s.p$group= paste(s.p$year, s.p$sex, s.p$timeperiod, sep="_")

#make year factor
s.p$year<- as.factor(s.p$year)

plot.sp= ggplot(s.p, aes(x=photo, y=surv, col=year, shape=sex, group=group))+
  geom_point() +geom_line()+  #geom_smooth(method="lm")+
  theme_classic(base_size = 18)+
  ylab("Survival proportion") +xlab("Photoperiod")

#sharpie manipulation
s.s= read.csv("Survival_Kingsolver1996_Fig1-2.csv")
s.s$group= paste(s.s$year, s.s$sex, s.s$timeperiod, sep="_")

#make year factor
s.s$year<- factor(s.s$year, levels=c("1991", "1992", "1993"), ordered=TRUE)

#omit longer final time period in both
s.s= s.s[-which(s.s$year==1991 &  s.s$timeperiod==3),]
s.s= s.s[-which(s.s$year==1993 &  s.s$timeperiod==4),]

plot.ss= ggplot(s.s, aes(x=treat, y=surv, col=year, shape=sex, group=group, factor(year)))+
  geom_point() +geom_line()+  #geom_smooth(method="lm")+
  theme_classic(base_size = 18)+
  ylab("Survival proportion") +xlab("treatment")

#Combine Figure 7
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/figures/")
pdf("Fig7_trait_behav.pdf",height = 5, width = 9)

plot.wt + plot.beh + 
  plot_layout(widths = c(1.5, 1.5))+
  plot_annotation(tag_levels = 'a')

dev.off()

#Figure 8
#combine selection data
names(s.p)[3]<-"treat"
s.s$experiment="Direct manipulation"
s.p$experiment="Photoperiod"
sel= rbind(s.s, s.p)
#change treatment terms
sel$pheno= "lighter"
sel$pheno[sel$treat %in% c("black","SD")]<-"darker"
sel$pheno= factor(sel$pheno, levels=c("lighter","darker"), ordered=TRUE)
#change sex labels
sel$sex[sel$sex=="both"]<-"combined"
#order year
sel$year<- factor(sel$year, levels=c("1991", "1992", "1993"), ordered=TRUE)
#order experiment
sel$experiment= factor(sel$experiment, levels=c("Photoperiod", "Direct manipulation"), ordered=TRUE)

plot.sel= ggplot(sel, aes(x=pheno, y=surv, col=year, shape=sex, group=group, factor(year)))+
  geom_point(size=3.5) +geom_line(lwd=1.5)+  #geom_smooth(method="lm")+
  facet_wrap(~experiment)+
  theme_classic(base_size = 20)+
  ylab("Survival proportion") +xlab("Wing phenotype")+
  scale_color_manual(values=c("#440154FF", "#238A8DFF","#B8DE29FF") )+
  scale_x_discrete(expand=c(0.1,0.1))

setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/figures/")
pdf("Fig8_selection.pdf",height = 6, width = 8)
plot.sel
dev.off()
#===================
#Figure 8. Plasticity and seasonality figure
#Extract photoperiod reaction norms from mean discriminant functions?
