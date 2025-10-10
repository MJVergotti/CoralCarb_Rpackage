#' Temperature conversion to kelvin
#'
#' Converts temperature in Celsius to kelvin
#' @param TCel The temperature in degrees Celsius
#' @return The temperature in kelvin
#' @examples
#' temp_kelvin <- TCel_to_Kelv(0);
#' temp_kelvin <- TCel_to_Kelv(c(0, -273));
#' @export
TCel_to_Kelv <- function(TCel){
  TKelv <- TCel + 273;
  return(TKelv);
}

#'Ksp <- 10^(-171.945 - 0.077993*TKelv + (2903.293/TKelv) + 71.595*log10(TKelv) + 
#'             ( - 0.068393 + 0.0017276*TKelv + 88.135/TKelv) * (Sal^0.5) - 0.10018*Sal + 0.0059415*Sal^1.5) # with your Salinity data in UPS.
#'
#'Btot <- 4.16*10^-4 * (Sal/35)*1000 # with your Salinity data in UPS. 
#'
#'K1 <- 10^-(3670.7/TKelv - 62.008 + 9.7944*log(TKelv) - 0.0118*Sal + 0.000116*Sal^2) # with your Salinity data in UPS.
#'
#'2 <- 10^-(1394.7/TKelv + 4.777 - 0.0184*Sal + 0.000118*Sal^2) # with your Salinity data in UPS.
#'
#'Kb <- 10^-(-log10(exp(((- 8966.9 - 2890.53*Sal^(1/2) - 77.942*Sal + 
#'                          1.728*Sal^(3/2) - 0.0996*Sal^(2))/(TCel + 273.15)) +
#'                        148.0248 + 137.1942*Sal^(1/2) + 1.62142*Sal - 
#'                        (24.4344 + 25.085*Sal^(1/2) + 0.2474*Sal)*log(TCel + 273.15) + 
#'                        0.053105*Sal^(1/2)*(TCel + 273.15))))
#'
#'pHcf <- - log10(Kb) - log10((d11Bsw - d11B)/(alpha_klochko*d11B - d11Bsw + 1000*(alpha_klochko - 1)))
#'
#'pHsw_cladocora <- (1/0.48)*pHcf - 10 
#'
#'protons <- 10^-(pHcf)
#'
#'Kd <- Kd_0 * exp( - K_kd*protons)*1000
#'
#'borate <- (Btot / (1 + (protons/Kb))) * 1000
#'
#'carbonate <- Kd*borate/(BCa/1000)


#' Dissolved inorganic carbon of the calcifying fluid
#'
#' Calculates the DIC of the calcifying fluid based on D'Olivo and McCulloch 2017
#' @param carbonate The concentration of carbonate ions in umol/kg. You can use the function carbonate() to calculate it. 
#' @param protons The concentration of hydrogen ions in mol/kg. You can use the function protons() to calculate it. 
#' @param K1 The first equilibrium constant of dissociation of carbonic acid in seawater. You can use the function K1() to calculate it. 
#' @param K2 The second equilibrium constant of dissociation of carbonic acid in seawater. You can use the function K2() to calculate it. 
#' @return The DICcf
#' @export
DICcf <- function(carbonate, protons, K1, K2){
  DICcf <- carbonate * (1 + (protons/K2) + (protons/K1*K2));
  return(DICcf);
}

#' Aragonite saturation state of the calcifying fluid
#'
#' Calculates the aragonite saturation state of the calcifying fluid using the formula from D'Olivo and McCUlloch 2017
#' @param carbonate The concentration of carbonate ions in umol/kg. You can use the function carbonate() to calculate it. 
#' @param ca_sw The concentration of Ca ions in mol/kg. Default is Ca_sw = 0.0112 (value for the Mediterranean Sea, based on McCulloch et al. 2012)
#' @param Ksp The solubility constant of aragonite in mol/kg. You can use the function Ksp() to calculate it. 
#' @return The aragonite saturation state of the calcifying fluid
#' @export
arg_sat_cf <- function(carbonate, ca_sw, Ksp){
  arg_sat_cf <- ((carbonate/(10^6))* ca_sw)/Ksp;
  return(arg_sat_cf);
}