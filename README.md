# CoralCarb
An R package to calculate the carbonate parameters of coral calcifying fluid. 

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R-CMD-check](https://github.com/MJVergotti/CoralCarb_Rpackage/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jdoe/CoralCarb/actions/workflows/R-CMD-check.yaml)

## Overview
`CoralCarb` provides a set of functions based on published equations in an attempt to facilitate coral researchers workflow when calculating carbonate chemistry parameters within the coral calcifying fluid, using measured biological inputs (&delta;<sup>11</sup>B and B/Ca) and environmental inputs (temperature and salinity). 

## Installation
The `CoralCarb` package is currently only available from GitHub. 
It can be installed by runing the following code:
```r
# If you don't have devtools installed:
# install.packages("devtools")

# Install the CoralCarb package
devtools::install_github("https://github.com/MJVergotti/CoralCarb_Rpackage.git")
```

## Use
The main functions are to calculate the parameters of the coral upregulation system (pH<sub>cf</sub>, DIC<sub>cf</sub>, &Omega;<sub>cf</sub>). 
It can also be used to apply the species specific calibration equations currently available in the literature to reconstruct pH<sub>sw</sub> from pH<sub>cf</sub>.

### Example: reconstruct pHsw from coral skeletons
```r
# Load the package
library(CoralCarb)

# 1. Calculate pHsw from the pHcf obtained from Cladocora caepsitosa using the following species-specific equation:
# Delta_pH_Cladocora = 4.7974 - 0.521*pHsw; with Delta_pH = pHcf - pHsw (Trotter et al. 2011)
estimated_pHcf <- 8.7
pHsw_Cladocora <- pHsw_coral(pHcf = estimated_pHcf, species = "Cladocora caespitosa")
print(pHsw_Cladocora) 

# 2. Calculate pHsw from the pHcf obtained from Porites sp. using the following species-specific equation:
# Delta_pH_Porites = 5.96 - 0.68*pHsw; with Delta_pH = pHcf - pHsw (Krief et al. 2010)
estimated_pHcf <- 8.7
pHsw_Porites <- pHsw_coral(pHcf = estimated_pHcf, species = "Porites sp.")
print(pHsw_Porites)
```

## Comments and support
Feel free to provide feedback!

* **Bugs/Issues:** Please report any bugs or suggest enhancements on the [Issue Tracker](https://github.com/MJVergotti/CoralCarb_Rpackage/issues).
* **Contact:** For questions, contact the package maintainer: <mj.vergotti@gmail.com>
