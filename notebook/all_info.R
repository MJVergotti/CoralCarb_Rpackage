# Convert temperature in degree celsius to temeprature in kelvin. 
TKelv <- TCel + 273 # with your temperature data in celsius.

# The solubility constant fo aragonite in mol/kg (Zeebe and Wolf-Gladrow 2001)
Ksp <- 10^(-171.945 - 0.077993*TKelv + (2903.293/TKelv) + 71.595*log10(TKelv) + 
             ( - 0.068393 + 0.0017276*TKelv + 88.135/TKelv) * (Sal^0.5) - 0.10018*Sal + 0.0059415*Sal^1.5) # with your Salinity data in UPS.

# Total boron concentration in mmol/kg (Zeebe and Wolf-Gladrow 2001)
Btot <- 4.16*10^-4 * (Sal/35)*1000 # with your Salinity data in UPS. 

# Equilibrium constants of dissotiation of carbonic acid in sw, using T in Kelvin 
# (Zeebe and Wolf-Gladrow 2001, based on Mehrback et al. 1973 and Dickson and Millero 1987)
K1 <- 10^-(3670.7/TKelv - 62.008 + 9.7944*log(TKelv) - 0.0118*Sal + 0.000116*Sal^2) # with your Salinity data in UPS.
K2 <- 10^-(1394.7/TKelv + 4.777 - 0.0184*Sal + 0.000118*Sal^2) # with your Salinity data in UPS.

# Constant of base dissociation constant
# (Dickson et al. 1999, Zeebe and Wolf-Gladrow 2001)
Kb <- 10^-(-log10(exp(((- 8966.9 - 2890.53*Sal^(1/2) - 77.942*Sal + 
                          1.728*Sal^(3/2) - 0.0996*Sal^(2))/(TCel + 273.15)) +
                        148.0248 + 137.1942*Sal^(1/2) + 1.62142*Sal - 
                        (24.4344 + 25.085*Sal^(1/2) + 0.2474*Sal)*log(TCel + 273.15) + 
                        0.053105*Sal^(1/2)*(TCel + 273.15)))) # with Sal in UPS and temperature in celsius. 


# the pH of the CF (based on Holcomb et al. 2014)
# You need to define a few constant first. 
# d11B, alpha 
# Here are some predefined:
# d11B in permil to choose from the following:
# d11Bsw_forster08 = 39.5 (default) (Foster 2008)
# d11Bsw_DOlivoMcCulloch17 = 39.61 (D'Olivo and McCUlloch 2017)
# alpha_klochko = 1.0272 (default) (Klochko et al. 2016)
pHcf <- - log10(Kb) - log10((d11Bsw - d11B)/(alpha_klochko*d11B - d11Bsw + 1000*(alpha_klochko - 1)))

# the pH of seawater
# intrude the a and b parameters of the calibration equation of your coral species
# Here the defalt is the Cladocora caespitosa equation defined by Trotter et al. 2011
# a_Trotter = 
# b_Trotter = 
pHsw_cladocora <- (1/0.48)*pHcf - 10 

# concentration of hydrogene ions in mol/kg
protons <- 10^-(pHcf))

# partition coefficient for B/Ca in mol/kg (McCulloch et al. 2017)
# define the constants:
# Kd_0 =
# K_kd = 
Kd <- Kd_0 * exp( - K_kd*protons)*1000

# concentration of borate ions in umol/kg (D'Olivo and McCulloch 2017, from Holcomb et al. 2014)
borate <- (Btot / (1 + (protons/Kb))) * 1000

# concentration of carbonate ions in umol/kg (based on McCulloch et al. 2017)
carbonate <- Kd*borate/(BCa/1000)

# amount of dissolved inorganic carbon in the CF (D'Olivo and McCulloch 2017)
DICcf <- carbonate * (1 + (protons/K2) + (protons/K1*K2))

# aragonite saturation state of the CF  (D'Olivo and McCulloch, 2012)
# define the constants:
# ca_sw_Med = 0.0112 # (default) value for the Mediterranean sea in mol/kg (McCulloch et al. 2012)
arg_sat_cf <- ((carbonate/(10^6))* ca_sw)/Ksp
