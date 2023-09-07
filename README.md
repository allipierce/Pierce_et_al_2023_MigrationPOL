
# Seasonal migration alters energetic trade-off optimization and shapes life history.

This repository contains the R model and analyses code for Pierce AP,
Yanco SW, & Wunder MB (2023) Seasonal migration alters energetic
trade-off optimization and shapes life history.

## Repository contents

### Model source code

Source code for the IBM and optimization model is available in the
[R](/R/) folder as a series of R scripts containing functions and
routines for model components as follows:

- [GAEMM.R](/GAEMM.R/)
- Top level function to run the metabolism and movement IBM and genetic
  algorithm optimization
- [expendenergy.R](/expendenergy.R/)
- R script containing DEB functions for the IBM metabolism sub-model
- [move.R](/move.R/)
- R script containing functions for the IBM movement sub-model
- [seasonalworld.R](/seasonalworld.R/)
- R script containing functions for simulating seasonal energy
  availability in the IBM