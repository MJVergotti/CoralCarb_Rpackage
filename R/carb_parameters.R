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


#' Aragonite solubility constant
#' 
#' Calculates the aragonite solubility constant in mol/kg based Zeebe and Wolf-Gladrow 2001
#' @param TKelv The temperature in Kelvin. You can use the function TCel() to calculate it 
#' @param Sal The salinity in UPS 
#' @return The aragonite solubility constant
#' @export
Ksp_arag <- function(TKelv, Sal){
  Ksp <- 10^(-171.945 - 0.077993*TKelv + (2903.293/TKelv) + 71.595*log10(TKelv) +
               ( - 0.068393 + 0.0017276*TKelv + 88.135/TKelv) * (Sal^0.5) - 0.10018*Sal + 0.0059415*Sal^1.5);
  return(Ksp);
}

#' Total boron concentration in the calcifying fluid
#'
#' Calculates the concentration of toal boron in the calcifying fluid in mmol/kg based on Zeebe and Wolf-Gladrow 2001
#' @param Sal The salinity in UPS
#' @return The total boron concentration
#' @export
Btot_cf <- function(Sal){
  Btot <- 4.16*10^-4 * (Sal/35)*1000;
  return(Btot);
}

#' 1st equilibrium constant of dissociation of carbonic acid in seawater 
#'
#' Calculates the 1st constant of dissociation of carbonic acid in seawater based on Zeebe and Wolf-Gladrow 2001, Mehrback et al. 1973 and Dickson and Millero 1987
#' @param TKelv The temperature in Kelvin. You can use the function TCel() to calculate it 
#' @param Sal The salinity in UPS
#' @return The 1st equilibrium constant of dissociation
#' @export
K1_sw <- function(TKelv, Sal){
  K1 <- 10^-(3670.7/TKelv - 62.008 + 9.7944*log(TKelv) - 0.0118*Sal + 0.000116*Sal^2);
  return(K1);
}

#' 2nd equilibrium constant of dissociation of carbonic acid in seawater
#' 
#' Calculated the 2nd constant of dissociation of carbonic acid in seawater based on Zeebe and Wolf-Gladrow 2001, Mehrback et al. 1973 and Dickson and Millero 1987
#' @param TKelv The temperature in Kelvin. You can use the function TCel() to calculate it 
#' @param Sal The salinity in UPS
#' @return The 2nd equilibrium constant of dissociation
#' @export
K2_sw <- function(TKelv, Sal){
  K2 <- 10^-(1394.7/TKelv + 4.777 - 0.0184*Sal + 0.000118*Sal^2);
  return(K2);
}

#' Constant of base dissociation of boron in seawater
#' 
#' Calculates the constant of base dissociation of boron in seawater based on Dickson et al. 1999, Zeebe and Wolf-Gladrow 2001
#' @param Sal The salinity in UPS
#' @param TCel The temperature in degrees Celsius
#' @return The constant of base dissociation
#' @export
Kb_sw <- function(Sal, TCel){
  Kb <- 10^-(-log10(exp(((- 8966.9 - 2890.53*Sal^(1/2) - 77.942*Sal + 
                            1.728*Sal^(3/2) - 0.0996*Sal^(2))/(TCel + 273.15)) +
                          148.0248 + 137.1942*Sal^(1/2) + 1.62142*Sal - 
                          (24.4344 + 25.085*Sal^(1/2) + 0.2474*Sal)*log(TCel + 273.15) + 
                          0.053105*Sal^(1/2)*(TCel + 273.15))));
  return(Kb);
}

#' Calcifying fluid pH
#' 
#' Calculates the pHcf based on Holcomb et al. 2014
#' @param Kb The constant of base dissociation. You can use the function Kb() to calculate it 
#' @param d11Bsw The d11B of seawater in permil. Default is d11B_foster08 = 39.5 (Foster 2008); else d11B_DOlivoMcCulloch17 = 39.61 (D'Olivo and McCUlloch 2017)
#' @param d11B The d11B of the coral
#' @param alpha_klochko The alpha value reported by Klochko: alpha = 1.0272 (Klochko et al. 2016)
#' @return The pHcf
#' @export
pHcf <- function(Kb, d11Bsw, d11B, alpha_klochko){
  pHcf <- - log10(Kb) - log10((d11Bsw - d11B)/(alpha_klochko*d11B - d11Bsw + 1000*(alpha_klochko - 1)));
  return(pHcf);
}

#' Seawater pH reconstruction
#' 
#' Calculates the reconstructed pH of seawater using the species specific equation for the coral Cladocora caespitosa, defined by Trotter et al. 2011. 
#' @param pHcf The pH of the calcifying fluid. You can use the function pHcf() to calculate it. 
#' @param species Choose a coral species from the list: 
#' "Porites cylindrica" (Hönisch et al. 2004), 
#' "Acropora nobilis" (Hönisch et al. 2004),
#' "Acropora sp." (Reynaud et al. 2004),
#' "Porites sp." (Krief et al. 2010),
#' "Stylophora pistillata" (Krief et al. 2010),
#' "Cladocora caespitosa" (Trotter al. 2011)
#' @return The pHsw
#' @export
pHsw_coral <- function(pHcf, species = "Cladocora caespitosa") {
  
  # Convert input variables to numeric
  pHcf <- as.numeric(pHcf)
  
  # Verify that the conversion was succesfull
  if(is.na(pHcf)){
    stop("Error: The argument 'pHcf' must be numeric or convertible to numeric.")
  }
  
  # Initialize variables
  c_val <- NULL
  d_val <- NULL
  
  # --- 1. Parameter Selection based on Species ---
  if (species == "Cladocora caespitosa") {
    # pHsw <- (pHcf - 4.7974) / (1 - 0.521)
    c_val <- - 0.521
    d_val <- 4.7974
    
  } else if (species == "Porites sp.") {
    # pHsw <- (pHcf - 5.96) / (1 - 0.68)
    c_val <- - 0.68
    d_val <- 5.96
    
  } else if (species == "Porites cylindrica") {
    # pHsw <- (pHcf - 4.72) / (1 - 0.53)
    c_val <- - 0.53
    d_val <- 4.72
    
    
  } else if (species == "Acropora nobilis") {
    # pHsw <- (pHcf - 4.40) / (1 - 0.50)
    c_val <- - 0.50
    d_val <- 4.40    
    
    
  } else if (species == "Acropora sp.") {
    # pHsw <- (pHcf - 4.28) / (1 - 0.48)
    c_val <- - 0.48
    d_val <- 4.28    
    
    
  } else if (species == "Stylophora pistillata") {
    # pHsw <- (pHcf - 6.06) / (1 - 0.70)
    c_val <- - 0.70
    d_val <- 6.06    
    
  } else {
    stop("Error: Species '", species, "' not recognized. Please use 'Cladocora caespitosa', 'Porites sp.', 'Porites cylindrica',
         'Acropora nobilis', 'Acropora sp.' or 'Stylophora pistillata'.")
  }
  # --- 2. Generic formula ---
  pHsw <- (pHcf - d_val) / (1 + c_val)
  
  return(pHsw)
}

#' Hydrogen ions concentrations in the calcifying fluid
#' 
#' Calculates the concentration of hydrogen ions in calcifying fluid in mol/kg
#' @param pHcf The pHcf. You can use the function pHcf() to calculate it 
#' @return The Hydrogen ions concentration
#' @export
protons_cf <- function(pHcf){
  protons <- 10^-(pHcf);
  return(pHcf);
}

#' Partition coefficient of B/Ca 
#' 
#' Calculates the partition coefficient of B/Ca in mol/kg (McCulloch et al. 2017)
#' @return The partition coefficient of B/Ca
#' @export
Kd_BCa <- function(protons){
  Kd <- 0.00297 * exp( - (- 0.0202)*protons)*1000;
  return(Kd);
}

#' Borate ion concentration in the calcifying fluid
#' 
#' Calculates the borate ion concentration in the calcifying fluid in umol/kg based on D'Olivo and McCulloch, 2017, Holcomb et al. 2014
#' @param Btot The total boron concentration. You can use the function Btot() to calculate it 
#' @param protons The concentration of hydrogen ions. You can use the function protons() to calculate it 
#' @param Kb The constant of base dissociation. You can use the function Kb() to calculate it 
#' @return The concentration of borate ions
#' @export
borate_cf <- function(Btot, protons, Kb){
  borate <- (Btot / (1 + (protons/Kb))) * 1000;
  return(borate);
}

#' Carbonate ion concentration in the calcifying fluid
#' 
#' Calculates the concentration of carbonate ions in the calcifying fluid in umol/kg (based on McCulloch et al. 2017)
#' @param Kd The partition coefficient of B/Ca. You can calculate it with the function carbonate()
#' @param borate The borate ion concentration. You can calculate it with the function borate()
#' @param BCa The B/Ca availability in the coral skeleton in mmol/mol
#' @return The carbonate ion concentration
#' @export
carbonate_cf <- function(Kd, borate, BCa){
  carbonate <- Kd*borate/(BCa/1000);
  return(carbonate);
}

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