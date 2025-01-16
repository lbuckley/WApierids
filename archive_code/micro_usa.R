# Adapt Colias niche model

#global historical climatology network?

#Seattle
#Center for urban horticulture
lat=47.657628
lon=-122.290255

dstart <- paste0("01/05/2017") # start date
dfinish <- paste0("01/09/2017") # end date

#Corfu site in central Washington in the Columbia National Wildlife Refuge
#12 miles west of Othello
#46.850663264 -119.535331192
  
#------------------

lonlat <- c(lon, lat) # (longitude, latitude)

DEP <- c(0, 3, 5, 10, 15, 20, 30, 50, 100, 200)

#----
library(NicheMapR)

micro <- micro_usa(loc = lonlat, dstart = dstart, dfinish = dfinish, DEP = DEP,
                   runmoist = 0, runshade = 0, Usrhyt = 0.01)

variable <- micro$metout[, "TALOC"]

variableDOY <- micro$metout[, "DOY"]
variableHOUR <- micro$metout[, "TIME"]

vals <- c()
begin <- 1 
end <- begin + 24*31
vals <- variable[begin:end]
valsDOY <- variableDOY[begin:end]
valsHOUR <- variableHOUR[begin:end]

days <- c()
for (i in 1:31) {
  days <- c(days, paste0("2017-0", month, "-", i))
}

df <- data.frame("Date" = rep(days, each = 24)[1 : (24 * 31)],
                 "Hour" = rep(0 : 23, 31)[1 : (24 * 31)],
                 "Data" = vals[0 : (24 * 31)])


df$Date <- format(as.POSIXct(paste0(df$Date, " ", df$Hour, ":00")), format = "%Y-%m-%d %H:%M")