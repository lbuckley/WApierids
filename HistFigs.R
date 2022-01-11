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

#---------
#Kingsolver and Gomulkiewicz. 2003. Environmental variation and selection. https://doi.org/10.1093/icb/43.3.470
#Consumption rate 24h

#Figure 1
#6hr growth rate for 4th and 5th instar
#4th and 5th instar
temps= c(11,17,23,29,35,40)
rgr4= c(0.004, 0.0072, 0.0195, 0.0224, 0.0325, 0.0174)

temps= c(11,17,23,29,35,41)
rgr5= c(0.00535, 0.011, 0.0176, 0.022, 0.0263, 0.0015)

#---------------------
#density plot
p1= ggplot(to, aes(x=Tmodel))+
  geom_density(lwd=1.2)+
  xlim(0,45)+
  xlab("Temperature (°C)")+
  theme_bw(base_size = 18)+theme(legend.position = c(0.6, 0.7))

# add selection arrows
#i + geom_segment(aes(x = 5, y = 30, xend = 3.5, yend = 25),
#                 arrow = arrow(length = unit(0.5, "cm")))
p2= p1 + geom_segment(data=sg, aes(x = temps, y = .09, xend = temps, yend = .09+pm.sg/20),
                      arrow = arrow(length = unit(0.5, "cm")))

#add grwoth rate
p3= p2 + geom_line(data=gr, aes(x=temps, y= mgr), color="blue", lwd=1.2)

#add consumption
p4= p3 + geom_line(data=cr, aes(x=temps, y= cr), color="green", lwd=1.2)

#----------------------
#plot out
setwd("/Volumes/GoogleDrive/My Drive/Buckley/Work/Proposals/NSF_ORCC/figures/")
pdf("TPCplot.pdf",height = 6, width = 10)

dev.off()

#===================
#Figure 7. Field selection figure



#===================
#Figure 8. Plasticity and seasonality figure
#Extract photoperiod reaction norms from mean discriminant functions?
