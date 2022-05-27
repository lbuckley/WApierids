#' @title Operative Environmental Temperature of a Butterfly
#'
#' @description Predicts body temperatures (operative environmental temperatures) of a butterfly in C.
#'  \cr \cr
#'  Based on \insertCite{Kingsolver1983;textual}{TrenchR} and \insertCite{Buckley2012;textual}{TrenchR}.
#' 
#' @param T_a \code{numeric} air temperature in C.
#'
#' @param Tg  \code{numeric} surface temperature in C in the sunlight.
#'
#' @param Tg_sh \code{numeric} surface temperature in C in the shade.
#'
#' @param u \code{numeric} wind speed in m / s.
#'
#' @param H_sdir \code{numeric} direct solar radiation flux in W / \ifelse{html}{\out{m<sup>2</sup>}}{\eqn{m^2}{ASCII}}.
#'
#' @param H_sdif \code{numeric} diffuse solar radiation flux in W / \ifelse{html}{\out{m<sup>2</sup>}}{\eqn{m^2}{ASCII}}.
#'
#' @param z \code{numeric} solar zenith angle in degrees.
#'
#' @param D \code{numeric} thoracic diameter in cm.
#'
#' @param delta \code{numeric} thoracic fur thickness in mm.
#'
#' @param alpha \code{numeric} wing solar absorbtivity as a proportion. Range for Colias is 0.4 to 0.7.
#'
#' @param r_g \code{numeric} substrate solar reflectivity (proportion), see \insertCite{Kingsolver1983;textual}{TrenchR}.
#' 
#' @param wing_angle \code{numeric} wing angle (degree)
#'
#' @param shade \code{logical} indicator whether body temperature should be calculate in sun (\code{FALSE}) or shade (\code{TRUE}).
#'
#' @return \code{numeric} predicted body (operative environmental) temperature (C).
#'
#' @details 
#'  Thermal radiative flux is calculated following \insertCite{Gates1980;textual}{TrenchR} based on \insertCite{Swinbank1960;textual}{TrenchR} and \insertCite{Kingsolver1983;textual}{TrenchR} estimates using the Brunt equation with black body sky temperature from \insertCite{Swinbank1963;textual}{TrenchR}.
#'
#' @family biophysical models
#'
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
#' @examples
#'   Tb_butterfly(
#'     T_a = 25, 
#'     Tg = 25, 
#'     Tg_sh = 20, 
#'     u = 0.4, 
#'     H_sdir = 300, 
#'     H_sdif = 100, 
#'     z = 30, 
#'     D = 0.36, 
#'     delta = 1.46, 
#'     alpha = 0.6, 
#'     r_g = 0.3)
#'
Tb_butterfly <- function (T_a, Tg, Tg_sh, u, H_sdir, H_sdif, z, D, delta, alpha, r_g = 0.3, wing_angle=42, shade = FALSE) {
  
  stopifnot(u >= 0, H_sdir >= 0, H_sdif >= 0, z >= -90, z <= 90, D > 0, delta >= 0, alpha >= 0, r_g >= 0, r_g <= 1, shade %in% c(FALSE, TRUE) )  
  
  # conversions
  
  # temperatures C to K
  
  TaK <- T_a + 273.15
  TaK_sh <- TaK
  Tg <- Tg + 273.15 
  Tg_sh <- Tg_sh + 273 
  
  # wind speed m/s to cm/s
  
  u <- u *100  
  
  # solar flux W/m2 to mW/cm2
  
  H_sdir <- H_sdir / 10 
  H_sdif <- H_sdif / 10
  
  # thoracic fur thickness mm to cm
  
  delta <- delta / 10     
  
  
  # Total solar radiation
  
  H_sttl <- H_sdir + H_sdif
  
  
  # Butterfly Parameters
  
  # surface emisivity, ranges from 0.95-1
  
  epsilon_s <- 0.97; 
  
  # Stefan-Boltzman constant, mW cm^-2 K^04 or 5.67*10^-8 W m-2 K-4
  
  sigma <- 5.67*10^-9 
  
  # butterfly thermal emmissivity
  
  Ep <- 1 
  
  # k_e- thermal conductivity of the fur, 1.3mWcm^-1*K^-1
  
  k_e <- 1.3; 
  
  # r_i- body radius from Kingsolver 1983
  
  r_i <- 0.15 
  
  # approximate thermal conductivity of air, mWcm^-1*K^-1
  
  k_a <- 0.25
  
  #  kinematic viscocity of air, cm^2/s at 300K 
  #     http://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
  
  v <- 15.68 * 10^-2  
  
  #-----
  #https://www.jstor.org/stable/pdf/4217669.pdf
  #S #area for radiation intercepted
  wing.angle= wing_angle/2*pi/180 #angle from vertical in radians, angle for pontia during flight, wing spread angle of 42degrees
  #radiation intercepted
  OC= (r+r/tan(wing.angle))*(1+tan(2*wing.angle)/tan(wing.angle))
  S= 2*OC*sin(wing.angle)
  S= S*alpha #account for reflectance
  
  L= 0.915 #calculate area for solar radiation using FWL, look for body length
  
  # convert length to area 
  # how to calculate absorptivity from wing traits
  #------
  
  A_sdir=S*L
  
  # total surface area cm^2 as area of cylinder without ends
  A_sttl <- pi * D * 2 
  
  #-----------------------------
  #UPDATE RADIATION ESTIMATES
  
  # direct and reflected surface areas For butterflies basking with wings perpendicular to radiation 
  A_sref <- 0 #CHECK A_sdir
  
  # RADIATIVE HEAT FLUx, mW
  
  #For Colias: Q_s <- alpha * A_sdir * H_sdir / cos(z * pi / 180) + alpha * A_sref * H_sdif + alpha * r_g * A_sref * H_sttl  
  
  Q_s <- A_sdir * H_sdir / cos(z * pi / 180) + alpha * A_sref * H_sdif + alpha * r_g * A_sref * H_sttl  
  #direct radiation already accounts for absorptivity
  
  #-----------------------------
  
  # THERMAL RADIATIVE FLUX in K
  
  # black body sky temperature from Swinbank 1963
  
  # Tsky <- 0.0552*(TaK)^1.5 
  
  # Gates 1980 Biophysical ecology based on Swinbank 1960, Kingsolver (1983) estimates using Brunt equation
  
  Tsky <- (1.22 * T_a - 20.4) + 273.15 
  
  # Q_t <- 0.5* A_sttl * Ep * sigma * (Tb^4 - Tsky^4) +0.5* A_sttl * Ep * sigma * (Tb^4 - Tg^4)
  
  
  # CONVECTIVE HEAT FLUX
  
  # Reynolds number- ratio of interval viscous forces
  
  R_e <- u * D / v
  
  # Nusselt number- dimensionless conductance
  
  N_u <- 0.6 * R_e^0.5
  
  # Kingsolver 1983
  
  # N_u <- 2.3 
  
  
  h_c <- N_u * k_a / D;
  
  
  # total convective heat transfer coefficient
  
  h_T <- (1 / h_c + (r_i + delta) * log((r_i + delta) / r_i) / k_e)^-1;  
  
  # convective heat transfer surface area
  
  # A_c <- A_sttl
  
  # Q_c <- h_T * A_c * (Tb-T_a)
  
  
  # Shade Adjustments
  
  ## UPDATE SHADE
  # if (shade) {
  #   
  #   # Calculate without basking by dividing areas by two
  #   
  #   A_sttl <- A_sttl / 2
  #   
  #   # RADIATIVE HEAT FLUX IN SHADE, mW
  #   
  #   A_sdir <- A_sttl/2
  #   A_sref <- A_sdir
  #   
  #   # No direct radiation, only diffuse and reflected
  #   
  #   H_sdir_sh <- 0
  #   H_sdif_sh <- H_sdif
  #   H_sttl <- H_sdif + H_sdif_sh 
  #   
  #   Q_s <- alpha * A_sdir * H_sdir_sh / cos(z * pi / 180) + alpha * A_sref * H_sdif_sh + alpha * r_g * A_sref * H_sttl; 
  #   
  #   # Use shaded surface temperature
  #   
  #   Tg< - Tg_sh
  #   
  # }
  
  # Solution 
  
  a <- A_sttl * Ep * sigma
  b <- h_T * A_sttl
  d <- h_T * A_sttl * TaK +0.5 * A_sttl * Ep * sigma * Tsky^4 + 0.5 * A_sttl * Ep * sigma * (Tg)^4 + Q_s
  
  # in K
  
  Te <- 1 / 2 * sqrt((2 * b) / (a * sqrt((sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1 / 3) / (2^(1 / 3) * 3^(2 / 3) * a) - (4 * (2 / 3)^(1 / 3) * d) / (sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1 / 3))) - (sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1 / 3) / (2^(1 / 3) * 3^(2 / 3) * a) + (4 * (2 / 3)^(1 / 3) * d) / (sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1 / 3)) - 1 / 2 * sqrt((sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1 / 3) / (2^(1 / 3) * 3^(2 / 3) * a) - (4 * (2 / 3)^(1 / 3) * d) / (sqrt(3) * sqrt(256 * a^3 * d^3 + 27 * a^2 * b^4) + 9 * a * b^2)^(1 / 3)) 
  
  # in C
  
  Te - 273.15
  
} 

#Plasticity: https://doi.org/10.1093/icb/38.3.545
#Kingsolver, J. G. 1987. Evolution and coadaptation of thermoregulatory behavior and wing pigmentation
#pattern in pierid butterflies. Evolution 41:472—490.
#P. occidentalis
#Longday= 16 hr light; shortday= 10 hr light


#----
#Thermoregulatory significance of wing melanization in Pieris butterflies: Physics, posture, and pattern
#https://www.jstor.org/stable/pdf/4217669.pdf
#reflectance model
#a dimensionless wing length, (L/D), the ratio of thermally-effective wing length to body diameter, and a dimensionless amount of radiation intercepted, (S/D). 

#Measurements with a spectroreflectometer (Kingsolver 1983 a) for Pieris show that the white dorsal wing surfaces have a solar reflectivity of more than 0.8, while the reflectivity of the black melanic wing regions is less than 0.3. 

#Thermal ecology of Pieris butterflies (Lepidoptera: Pieridae): a new mechanism of behavioral thermoregulation
#https://www.jstor.org/stable/pdf/4217668.pdf
#Colorado temperature dependence of flight

orient.angle=0 #orientation angle, angle from vertical

wing.angle.d= 45 #angle from vertical
#convert to radians
wing.angle= pi/10   #wing.angle.d*pi/180

#basking during flight initiation: Pontia (42degrees), Pieris (17 degrees)
#basking during flight behavior: Pontia (48degrees), Pieris (23 degrees)

#body radius
r= 0.36/2 #cm CHECK
#wingspan 4.5 to 6.5 cm #https://entnemdept.ufl.edu/creatures/veg/leaf/imported_cabbageworm.htm

R=0.8 #reflectivity
#L #reflective wing length
#S #area for radiation intercepted
#r #body radius

#reflective wing length
L= (r+r/tan(wing.angle))*tan(2*wing.angle)/tan(wing.angle)+r
#check degree vs radians

#radiation intercepted
OC= (r+r/tan(wing.angle))*(1+tan(2*wing.angle)/tan(wing.angle))
S= 2*OC*sin(wing.angle)
S=S*R #check reduction
#For wings of reflectivity R less than 1.0, the radiation intercepted must be reduced appropriately from the values
#of S obtained above. For each reflection from the wing, the radiation intercepted is reduced by a factor R

#S is wing area absorbing direct solar radiation

#Recreate figure 3
#Fig 3: maximum radiation intercepted (S/D), minimum wing length required to intercept maximum radiation (L/D) 

n=c(0, 1, 2, 3, 4, 5)
wing.a= pi/(2*(1+2*n))

Ss= rep(NA, length(wing.a))
Ls= rep(NA, length(wing.a))

for(wing.k in 1: length(wing.a)){
  wing.angle= wing.a[wing.k]
  
  #reflective wing length
  Ls[wing.k]= (r+r/tan(wing.angle))*tan(2*wing.angle)/tan(wing.angle)+r
  
  #radiation intercepted
  OC= (r+r/tan(wing.angle))*(1+tan(2*wing.angle)/tan(wing.angle))
  Ss[wing.k]= 2*OC*sin(wing.angle)
}

D=0.48 #to produce Fig. 3 #0.36
plot(wing.a*180/pi, Ss/D, type="l")
points(wing.a*180/pi, Ls/D, type="l", lty="dashed")

#plot lengths
plot(wing.a*180/pi, Ss, type="l")
points(wing.a*180/pi, Ls, type="l", lty="dashed")

#-----
#Fig 5
ref= c(0,.2,.4,.6,.8)
plot(wing.a*180/pi, Ss/D, type="l", ylim=c(0,6))
for (ref.k in 1:length(ref)){
points(wing.a*180/pi, Ss*ref[ref.k]/D, type="l")
}

#-----
#Fig 7
wing.a= seq(5,80,2)
Tbs= rep(NA, length(wing.a))

for(wing.k in 1: length(wing.a)){
  Tbs[wing.k]= Tb_butterfly( T_a = 20, Tg = 20, Tg_sh = 20, u = 1, 
                            H_sdir = 800, H_sdif = 200, z = 30, D = 0.36, 
                            delta = 1.46, alpha = 0.6, r_g = 0.3, wing_angle=wing.a[wing.k])
}
  
plot(wing.a,Tbs-20, xlab="Wing angle (degrees)", type="l")
  
