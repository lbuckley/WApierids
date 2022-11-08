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
dat.sub= subset(dat.all, dat.all$SOLR>0)

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
p1= ggplot(dat.sub, aes(x=TALOC))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Temperature at plant height (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 16)+theme(legend.position = c(0.45, 0.7))+
  geom_vline(xintercept=34.5)+
  ylim(0,0.072)+
  xlim(0,55)

p1.ref= ggplot(dat.sub, aes(x=TAREF))+
  geom_density(alpha=0.4, aes(fill=period, color=period))+  
  facet_wrap(~seas)+
  xlab("Temperature at reference height (°C)")+
  ylab("Density" )+
  theme_classic(base_size = 16)+theme(legend.position = c(0.45, 0.7))+
  geom_vline(xintercept=34.5)+
  ylim(0,0.072)+
  xlim(0,55)

#----------------
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/PlastEvolAmNat/figures/")
pdf("Fig_Pontia.pdf",height = 14, width = 6)
plot.sv2 / p1.ref / p1 + plot_layout(widths = c(1, 1.6))+ 
  plot_annotation(tag_levels = 'A')
dev.off()

