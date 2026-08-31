#' Temperature conversion to Kelvin
#'
#' Converts temperature in Celsius to Kelvin.
#' @param TCel The temperature in degrees Celsius.
#' @return The temperature in Kelvin.
#' @examples
#' temp_kelvin1 <- TCel_to_Kelv(0)
#' temp_kelvin2 <- TCel_to_Kelv(c(0, -273.15))
#' @export
TCel_to_Kelv <- function(TCel){
  TKelv <- TCel + 273.15
  return(TKelv)
}


#' Aragonite solubility constant
#' 
#' Calculates the aragonite solubility constant in mol/kg based Zeebe and Wolf-Gladrow 2001.
#' @param TCel The temperature in degrees Celsius.
#' @param Sal The salinity in UPS.
#' @return The aragonite solubility constant.
#' @examples
#' Ksp1 <- Ksp_arag(24, 38)
#' Ksp2 <- Ksp_arag(c(23, 25, 27), c(37, 37.5, 38))
#' @export
Ksp_arag <- function(TCel, Sal){
  Ksp <- 10^(-171.945 - 0.077993*(TCel + 273.15) + (2903.293/(TCel + 273.15)) + 71.595*log10((TCel + 273.15)) +
               ( - 0.068393 + 0.0017276*(TCel + 273.15) + 88.135/(TCel + 273.15)) * (Sal^0.5) - 0.10018*Sal + 0.0059415*Sal^1.5);
  return(Ksp)
}

#' Total boron concentration
#'
#' Calculates the total boron concentration in mmol/kg based on Zeebe and Wolf-Gladrow 2001.
#' @param Sal The salinity in UPS.
#' @return The total boron concentration.
#' @examples
#' Btot1 <- Btot(35)
#' Btot2 <- Btot(c(35, 38))
#' @export
Btot <- function(Sal){
  Btot <- 4.16*10^-4 * (Sal/35)*1000
  return(Btot)
}

#' 1st equilibrium constant of dissociation of carbonic acid in seawater 
#'
#' Calculates the 1st constant of dissociation of carbonic acid in seawater based on Zeebe and Wolf-Gladrow 2001, Mehrback et al. 1973 and Dickson and Millero 1987.
#' @param TCel The temperature in degrees Celsius.
#' @param Sal The salinity in UPS.
#' @return The 1st equilibrium constant of dissociation.
#' @examples
#' K1_1 <- K1_sw(25, 35)
#' K1_2 <- K1_sw(c(27, 28, 29), c(38.7, 38.6, 38.2))
#' @export
K1_sw <- function(TCel, Sal){
  K1 <- 10^-(3670.7/(TCel + 273.15) - 62.008 + 9.7944*log((TCel + 273.15)) - 0.0118*Sal + 0.000116*Sal^2)
  return(K1)
}

#' 2nd equilibrium constant of dissociation of carbonic acid in seawater
#' 
#' Calculated the 2nd constant of dissociation of carbonic acid in seawater based on Zeebe and Wolf-Gladrow 2001, Mehrback et al. 1973 and Dickson and Millero 1987.
#' @param TCel The temperature in degrees in Celsius.
#' @param Sal The salinity in UPS.
#' @return The 2nd equilibrium constant of dissociation.
#' @examples
#' K2_1 <- K2_sw(25, 35)
#' K2_2 <- K2_sw(c(27, 28, 29), c(38.7, 38.6, 38.2))
#' @export
K2_sw <- function(TCel, Sal){
  K2 <- 10^-(1394.7/(TCel + 273.15) + 4.777 - 0.0184*Sal + 0.000118*Sal^2)
  return(K2)
}

#' Constant of base dissociation of boron in seawater
#' 
#' Calculates the constant of base dissociation of boron in seawater based on Dickson et al. 1999, Zeebe and Wolf-Gladrow 2001.
#' @param Sal The salinity in UPS.
#' @param TCel The temperature in degrees Celsius.
#' @return The constant of base dissociation.
#' @examples
#' Kb_1 <- Kb_sw(25, 35)
#' Kb_2 <- Kb_sw(c(27, 28, 29), c(38.7, 38.6, 38.2))
#' @export
Kb_sw <- function(TCel, Sal){
  Kb <- 10^-(-log10(exp(((- 8966.9 - 2890.53*Sal^(1/2) - 77.942*Sal + 
                            1.728*Sal^(3/2) - 0.0996*Sal^(2))/(TCel + 273.15)) +
                          148.0248 + 137.1942*Sal^(1/2) + 1.62142*Sal - 
                          (24.4344 + 25.085*Sal^(1/2) + 0.2474*Sal)*log(TCel + 273.15) + 
                          0.053105*Sal^(1/2)*(TCel + 273.15))))
  return(Kb)
}

#' Calcifying fluid pH
#' 
#' Calculates the pHcf based on Holcomb et al. 2014.
#' @param Kb The constant of base dissociation. You can use the function Kb() to calculate it.
#' @param d11Bsw The d11B of seawater in permil. Default is d11Bsw = 39.5 (Foster 2008); else d11Bsw = 39.61 (D'Olivo and McCUlloch 2017).
#' @param d11B The d11B of the coral.
#' @param alpha_klochko The alpha value reported by Klochko: alpha = 1.0272 (Klochko et al. 2016).
#' @return The pHcf.
#' @examples
#' pHcf <- pHcf(2.6e-09, 39.5, 24.43, 1.0272)
#' @export
pHcf <- function(Kb, d11B, d11Bsw = 39.5, alpha_klochko = 1.0272){
  pHcf <- - log10(Kb) - log10((d11Bsw - d11B)/(alpha_klochko*d11B - d11Bsw + 1000*(alpha_klochko - 1)))
  return(pHcf)
}

#' Seawater pH reconstruction
#' 
#' Calculates the reconstructed pH of seawater using species-specific equations.
#' @param pHcf The pH of the calcifying fluid. You can use the function pHcf() to calculate it.
#' @param species Choose a coral species from the list: 
#' "Porites cylindrica" (Hönisch et al. 2004), 
#' "Acropora nobilis" (Hönisch et al. 2004),
#' "Acropora sp." (Reynaud et al. 2004),
#' "Porites sp." (Krief et al. 2010),
#' "Stylophora pistillata" (Krief et al. 2010),
#' "Cladocora caespitosa" (Trotter et al. 2011) (Default).
#' @return The pHsw.
#' @examples
#' pHsw_Cladocora <- pHsw_coral(8.51, species = "Cladocora caespitosa")
#' pHsw_default <- pHsw_coral(8.51)
#' @export
pHsw_coral <- function(pHcf, 
                       species = c("Cladocora caespitosa", 
                                   "Porites sp.", 
                                   "Porites cylindrica", 
                                   "Acropora nobilis", 
                                   "Acropora sp.", 
                                   "Stylophora pistillata")) {
  
  # Automatically validates user input against the allowed options above.
  # Defaults to the first item ("Cladocora caespitosa") if omitted.
  species <- match.arg(species)
  
  # Assign parameters based on species
  if (species == "Cladocora caespitosa") {
    c_val <- -0.521
    d_val <- 4.7974
  } else if (species == "Porites sp.") {
    c_val <- -0.68
    d_val <- 5.96
  } else if (species == "Porites cylindrica") {
    c_val <- -0.53
    d_val <- 4.72
  } else if (species == "Acropora nobilis") {
    c_val <- -0.50
    d_val <- 4.40
  } else if (species == "Acropora sp.") {
    c_val <- -0.48
    d_val <- 4.28
  } else if (species == "Stylophora pistillata") {
    c_val <- -0.70
    d_val <- 6.06
  } 
  
  # --- 2. Generic formula ---
  pHsw <- (pHcf - d_val) / (1 + c_val)
  
  return(pHsw)
}

#' Hydrogen ions concentrations in the calcifying fluid
#' 
#' Calculates the concentration of hydrogen ions in calcifying fluid in mol/kg.
#' @param pHcf The pHcf. You can use the function pHcf() to calculate it. 
#' @return The Hydrogen ions concentration.
#' @examples
#' H_1 <- protons_cf(8.51)
#' @export
protons_cf <- function(pHcf){
  protons <- 10^-(pHcf)
  return(protons)
}

#' Partition coefficient of B/Ca 
#' 
#' Calculates the partition coefficient of B/Ca in mol/kg (McCulloch et al. 2017).
#' @param protons The concentration of hydrogen ions. Calculate with the function protons().
#' @return The partition coefficient of B/Ca.
#' @examples
#' Kd_1 <- Kd_BCa(3.04e-09)
#' @export
Kd_BCa <- function(protons){
  Kd <- 0.00297 * exp( - (- 0.0202)*protons)*1000
  return(Kd)
}

#' Borate ion concentration in the calcifying fluid
#' 
#' Calculates the borate ion concentration in the calcifying fluid in umol/kg based on D'Olivo and McCulloch, 2017, Holcomb et al. 2014.
#' @param Btot The total boron concentration. You can use the function Btot() to calculate it.
#' @param protons The concentration of hydrogen ions. You can use the function protons() to calculate it.
#' @param Kb The constant of base dissociation. You can use the function Kb() to calculate it.
#' @return The concentration of borate ions.
#' @examples
#' BOH3_cf <- borate_cf(0.45, 3.04e-09, 2.58e-09)
#' @export
borate_cf <- function(Btot, protons, Kb){
  borate <- (Btot / (1 + (protons/Kb))) * 1000
  return(borate)
}

#' Carbonate ion concentration in the calcifying fluid
#' 
#' Calculates the concentration of carbonate ions in the calcifying fluid in umol/kg (based on McCulloch et al. 2017).
#' @param Kd The partition coefficient of B/Ca. You can calculate it with the function carbonate().
#' @param borate The borate ion concentration. You can calculate it with the function borate().
#' @param BCa The B/Ca availability in the coral skeleton in mmol/mol.
#' @return The carbonate ion concentration.
#' @examples
#' CO32_cf <- carbonate_cf(2.97, 208.3, 664.8)
#' @export
carbonate_cf <- function(Kd, borate, BCa){
  carbonate <- Kd*borate/(BCa/1000)
  return(carbonate)
}

#' Dissolved inorganic carbon of the calcifying fluid
#'
#' Calculates the DIC of the calcifying fluid based on D'Olivo and McCulloch 2017.
#' @param carbonate The concentration of carbonate ions in umol/kg. You can use the function carbonate() to calculate it. 
#' @param protons The concentration of hydrogen ions in mol/kg. You can use the function protons() to calculate it. 
#' @param K1 The first equilibrium constant of dissociation of carbonic acid in seawater. You can use the function K1() to calculate it. 
#' @param K2 The second equilibrium constant of dissociation of carbonic acid in seawater. You can use the function K2() to calculate it. 
#' @return The DICcf
#' @examples
#' DICcf_1 <- DICcf(930.5, 3.04e-09, 1.47e-06, 1.16e-09)
#' @export
DICcf <- function(carbonate, protons, K1, K2){
  DICcf <- carbonate * (1 + (protons/K2) + (protons/K1*K2))
  return(DICcf)
}

#' Aragonite saturation state of the calcifying fluid
#'
#' Calculates the aragonite saturation state of the calcifying fluid using the formula from D'Olivo and McCUlloch 2017.
#' @param carbonate The concentration of carbonate ions in umol/kg. You can use the function carbonate() to calculate it. 
#' @param ca_sw The concentration of Ca ions in mol/kg. Default is Ca_sw = 0.0112 (value for the Mediterranean Sea, based on McCulloch et al. 2012).
#' @param Ksp The solubility constant of aragonite in mol/kg. You can use the function Ksp() to calculate it. 
#' @return The aragonite saturation state of the calcifying fluid.
#' @examples
#' omega_cf <- arg_sat_cf(930.5, 0.0112, 7.25e-07)
#' @export
arg_sat_cf <- function(carbonate, Ksp, ca_sw = 0.0112){
  arg_sat_cf <- ((carbonate/(10^6))* ca_sw)/Ksp
  return(arg_sat_cf)
}

#' Master function that calculated all parameters of Coral Calcifying Fluid and Seawater Chemistry
#' 
#' Calculates the complete set of CF carbonate chemistry and seawater pH reconstruction from parameters from raw experimental and environmental data in a single step.
#' 
#' @param TCel Temperature in degrees Celsius.
#' @param Sal Salinity in UPS.
#' @param d11B The d11B of the coral skeleton (in permil).
#' @param BCa The B/Ca ratio of the coral skeleton (in mmol/mol).
#' @param d11Bsw Seawater d11B (in permil). Default is 39.5 (Foster 2008). Alternative is 39.61 (D'Olivo and McCulloch 2017).
#' @param alpha_klochko Isotopic fractionation factor alpha. Default is 1.0272 (Klochko et al. 2006).
#' @param species Coral species for pH reconstruction. Defaults to "Cladocora caespitosa".
#' @param ca_sw Concentration of Ca ions in seawater (mol/kg). Default is 0.0112 (McCulloch et al. 2012).
#' 
#' @return A \code{data.frame} (or vector if single inputs are provided) containing all calculated 
#'   constants, pH values, ion concentrations, and saturation states.
#' @examples
#' # Calculate for single values
#' coralCF(TCel = 25, Sal = 35, d11B = 24.43, BCa = 664.8)
#' 
#' # Vectorized calculation for multiple samples
#' df_results <- coralCF(
#'   TCel = c(25, 26), 
#'   Sal = c(35, 36), 
#'   d11B = c(24.43, 25.10), 
#'   BCa = c(664.8, 670.2)
#' )
#' @export
coralCF <- function(TCel, 
                    Sal, 
                    d11B, 
                    BCa, 
                    d11Bsw = 39.5, 
                    alpha_klochko = 1.0272, 
                    species = c("Cladocora caespitosa",
                                "Porites sp.", 
                                "Porites cylindrica", 
                                "Acropora nobilis", 
                                "Acropora sp.", 
                                "Stylophora pistillata"), 
                    ca_sw = 0.0112) {
  # Validate species selection
  species <- match.arg(species)
  
  # 1. Temperature & Constants
  TKelv <- TCel_to_Kelv(TCel)
  Ksp   <- Ksp_arag(TCel, Sal)
  B_tot <- Btot(Sal)
  K1    <- K1_sw(TCel, Sal)
  K2    <- K2_sw(TCel, Sal)
  Kb    <- Kb_sw(TCel, Sal)
  
  # 2. Calcifying Fluid pH & Seawater pH Reconstruction
  pH_cf <- pHcf(Kb = Kb, d11B = d11B, d11Bsw = d11Bsw, alpha_klochko = alpha_klochko)
  pH_sw <- pHsw_coral(pHcf = pH_cf, species = species)
  
  # 3. Ions & Speciation
  H_cf       <- protons_cf(pH_cf)
  Kd         <- Kd_BCa(H_cf)
  BOH3       <- borate_cf(Btot = B_tot, protons = H_cf, Kb = Kb)
  CO3        <- carbonate_cf(Kd = Kd, borate = BOH3, BCa = BCa)
  DIC_cf     <- DICcf(carbonate = CO3, protons = H_cf, K1 = K1, K2 = K2)
  Omega_arag <- arg_sat_cf(carbonate = CO3, Ksp = Ksp, ca_sw = ca_sw)
  
  # Return structured data frame
  results <- data.frame(
    TCel         = TCel,
    Sal          = Sal,
    TKelv        = TKelv,
    d11B         = d11B,
    BCa          = BCa,
    Kb           = Kb,
    Ksp          = Ksp,
    Btot         = B_tot,
    K1           = K1,
    K2           = K2,
    pHcf         = pH_cf,
    pHsw         = pH_sw,
    protons_cf   = H_cf,
    Kd_BCa       = Kd,
    borate_cf    = BOH3,
    carbonate_cf = CO3,
    DICcf        = DIC_cf,
    arg_sat_cf   = Omega_arag
  )
  
  return(results)
}
