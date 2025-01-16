# WApierids
Analyses associated with a resurvey of WA Pierid butterflies

# GENERAL INFORMATION

This README.txt file was updated on August 20, 2024 by Lauren Buckley

## A. Paper associated with this archive 
Citation: Buckley LB and Kingsolver JG. Functional resurveys and models reveal the interplay of plasticity and evolution of Pierid butterflies in response to recent climate change

Brief abstract: The extent of contemporary evolution will be an important determinant of the organismal and biodiversity consequences of climate change. Organisms can respond to climate change via tracking through space or time, phenotypic plasticity, or evolution. A key unknown shaping these interacting responses is whether plasticity facilitates or hinders evolution. We synthesize two resurvey projects for Pierid butterflies that evaluate the interplay of plasticity and evolution in responses to climate change. The temperature dependence of larval development and growth constitutes an important mechanistic link between phenotypes and the environment. Adult wing melanization, which influences the absorption of solar radiation to alter body temperature and flight activity, responds plastically to development environments. We examine shifts over recent decades in microclimate experienced by larval and adult Pierid butterflies. We assess the implications of the microclimate shifts for selection on the temperature sensitivity of larval feeding rate. We estimate selection for warmer thermal optima for feeding, particularly in response to spring warming. Selection estimates partially correspond to observed evolution, but we also detect some evolutionary surprises. The research highlights the importance of considering interacting organismal responses to climate change.

## B. Originators

Lauren B. Buckley, Department of Biology, University of Washington, Seattle, WA 98195-1800, USA
Joel G. Kingsolver, Department of Biology, University of North Carolina, Chapel Hill, NC 27599, USA

## C. Contact information
Lauren Buckley
Department of Biology, University of Washington, Seattle, WA 98195-1800, USA
lbuckley@uw.edu

## D. Dates of data collection
No new data are collected, but see references for past data utilized. 

## E. Geographic Location(s) of data collection
Montrose Valley, CO, N 38.62, W 108.02, 1633m
Sacramento Valley, CA, N 38.44, W121.86, 19m
Seattle, WA

## F. Funding Sources 
This research was supported by the US National Science Foundation (DEB-1120062 to L.B.B. and J.G.K., IOS- 2222089 to L.B.B., IOS- 2222090 to J.G.K).

# ACCESS INFORMATION

## 1. Licenses/restrictions placed on the data or code
CC0 1.0 Universal (CC0 1.0)
Public Domain Dedication

## 2. Data derived from other sources
data/era5_micro_shade and data/era5_micro_shade are from https://cds.climate.copernicus.eu/user/register and processed using micro_era5 as described in the micro_era5.R file.

data/Higgins from Higgins et al. 2014

data/KingsolverPrapae from Kingsolver 2024

data/NielsenColias from Nielsen and Kingsolver 2020

## 3. Recommended citation for this data/code archive
Buckley LB and Kingsolver JG. 2024. Functional resurveys and models reveal the interplay of plasticity and evolution of Pierid butterflies in response to recent climate change. https://github.com/lbuckley/WApierids. 

Data and code will be uploaded to Dryad upon acceptance.

# DATA & CODE FILE OVERVIEW

This data repository consist of 5 data folders, 2 code scripts, and this README document, with the following data and code filenames and variables

## Data files and variables

1. Higgins/
FR1974.csv
Historic data from Higgins et al. 2014

Temp: Temperature (C)
FR: Feeding rate
Pop: Population
Year: Sampling year

totalFT.csv
Contemporary data from Higgins et al. 2014

ID: ID
Trial: trial number
Fam: Family ID
Temp: Temperature (C)
lnTMG: natural log of Total Mass Gain, (final mass/initial mass)/time, 1/s	
pop: population

HigginsTPC.csv
Parameter fits from Higgins et al. 2014 for the following performance curve:
tpc= function(T, Fmax,To, row, sigma) Fmax*exp(-exp(row*(T-To)-6)-sigma*(T-To)^2)

species: Colias species
year: sample year
population: population collection location
Fmax: maximum feeding rate
Topt: optimal temperature
row: thermal sensitivity of feeding above Topt
sigma: thermal sensitivity of feeding below Topt

2. KingsolverPrapae/PrapaeUW.Seln2.1999.Combineddata.OPUS2021.csv
Recnum	record number	
Species	species	
Site	Study site	
Year	year	
Study	study season (Spring or Summer)	
ID	individual ID	
mother	mother (full-sib group) ID	
Jdate.Hatch	hatch date, Julian date of year	
Plant.ID	plant ID number	
Time.2nd	time to start of 2nd instar	d
Time.3rd	time to start of 3rd instar	d
Time.4th	time to start of 4th instar	d
Time.5th	time to start of 5th instar	d
Time.Pupa	time to  pupation	d
Time.Eclos	time to eclosion	d
Mass.2nd	mass at start of 2nd instar	mg
Mass.5th	mass at start of 5th instar	mg
Mass.Pupa	mass at pupation	mg
Mass.Eclos	mass at eclosion	mg
Surv.Pupa	survival to pupation? Y or N	
Surv.Eclos	survival to eclosion ? Y or N	
Sex	sex (M or F)	
Eggs	number of eggs laid	

3. NielsenColias/
DevelopmentTime_Table2_AllenSmith.csv
Data from Table 2 of Allen and Smith 1958

Temperature: Temperature (C)
DT_egg: development time as eggs (days)
DT_larvae: development time as larvae (days)
DT_pupae: development time as pupae (days)
DT_total: total development time (days)

Nielsen_conteptSummary.csv
Data summarized from Nielsen and Kingsolver 2020

Larval_Photoperiod: hours of light per day at which caterpillars were raised
means: mean of percent reflectance at 650 nm
sd: standard deviation of percent reflectance at 650 nm
n: number of individuals at each photoperiod
se: standard error of percent reflectance at 650 nm

Nielsen_oldData.csv is Hoffman_1973_extracted_data.csv from Nielsen and Kingsolver 2020
Data extracted from figures 3 and 4 in:
 Hoffmann RJ. 1973 Environmental Control of Seasonal Variation in the Butterfly Colias eurytheme. I. Adaptive Aspects of a Photoperiodic Response. Evolution 27, 387. (doi:10.2307/2407302)
Data was extracted using the Figure_Callibration package (http://www.astro.physik.uni-goettingen.de/~hessman/ImageJ/Figure_Calibration/) for imageJ (v1.52a).

Photoperiod: hours of light per day at which caterpillars were raised, same for both figures. Given preciesly by Hoffmann
Measured_Photo: Estimate value of photoperiod from the image (x-axis). Given separately preceding each measurement. Comparison to Photoperiod gives an idea of the accuracy of the data extraction.
Reflectance: Estimated mean percent reflectance at 650 nm from figure 3 (y-axis)
RefMinusSE: Estimated mean reflectance - standard error from figure 3 (y-axis, position of lower error bar)
RefMinusSE: Estimated mean reflectance + standard error from figure 3 (y-axis, position of upper error bar)
UpperLength: Estimated mean forewing length (mm) from figure 4 (y-axis)
ULenMinusSE: Estimated mean forewing length - standard error from figure 4 (y-axis, position of lower error bar)
ULenPlusSE: Estimated mean forewing length + standard error from figure 4 (y-axis, position of upper error bar)
LowerLength: Estimated mean hindwing length (mm) from figure 4 (y-axis)
LLenMinusSE: Estimated mean hindwing length - standard error from figure 4 (y-axis, position of lower error bar)
LLenPlusSE: Estimated mean hindwing length + standard error from figure 4 (y-axis, position of upper error bar)

4. era5_micro_shade and era5_micro_sun/
Files are output of the micro_era function in NicheMapR, assuming shade (minshade=90, maxshade=100) or sun (minshade=25, maxshade=30) as described in micro_era5.R. See https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html.

dates - vector of dates (POSIXct, UTC)
DOY - day-of-year
TIME - time of day (mins)
TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
RHLOC - relative humidity (%) at local height (specified by 'Usrhyt' variable)
RH - relative humidity (%) at reference height (specified by 'Refhyt', 2m default)
VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)
VREF - wind speed (m/s) at reference height (specified by 'Refhyt', 2m default)
SNOWMELT - snowmelt (mm)
POOLDEP - water pooling on surface (mm)
PCTWET - soil surface wetness (%)
ZEN - zenith angle of sun (degrees - 90 = below the horizon)
SOLR - solar radiation (W/m2) (unshaded, adjusted for slope, aspect and horizon angle)
TSKYC - sky radiant temperature (°C)
DEW - dew fall (mm / h)
FROST - frost (mm / h)
SNOWFALL - snow predicted to have fallen (cm)
SNOWDEP - predicted snow depth (cm)
SNOWDENS - snow density (g/cm3)
D0cm ... - soil temperature (°C) at each of the 10 specified depths

## Code scripts and workflow
1. Figs2_3_5_ColiasLarvae.R: code for producing figures 2, 3, and 5 and associated supplementary figures

2. Figs4_PierisLarvae_beta.R: code for producing figure 4 and associated supplementary figures

Previous code used during development is in the archive_code folder

# SOFTWARE VERSIONS
R version 4.1.0 (2021-05-18)

Packages:
ggplot2_3.4.4
reshape2_1.4.4
reshape_0.8.8
viridisLite_0.4.0
patchwork_1.2.0
TrenchR_1.1.1
tidyverse_1.3.2
rTPC_1.0.2
nls.multstart_1.2.0
dplyr_1.1.2

# REFERENCES
Allen, W., and R. Smith. 1958. Some factors influencing the efficiency of Apanteles medicaginis Muesebeck (Hymenoptera: Braconidae) as a parasite of the alfalfa caterpillar, Colias philodice eurytheme Boisduval. Hilgardia 28:1–42.

Higgins JK, MacLean HJ, Buckley LB, Kingsolver JG. Geographic differences and microevolutionary changes in thermal sensitivity of butterfly larvae in response to climate. Functional Ecology. 2014; 28:982-9. https://www.jstor.org/stable/24033587

Kingsolver, JG. Historical data for feeding, growth and life history in Pieris rapae [Dataset]. Dryad. 2024; https://doi.org/10.5061/dryad.nzs7h44zz

Nielsen ME, Kingsolver JG. Compensating for climate change–induced cue‐environment mismatches: evidence for contemporary evolution of a photoperiodic reaction norm in Colias butterflies. Ecology letters. 2020;23(7):1129-36. 
