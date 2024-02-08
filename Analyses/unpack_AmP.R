library(R.matlab)
library(tidyverse)
allStat <- readMat("./data/AmPdata/allStat.mat")

labs <- readMat("./data/AmPdata/allLabel.mat")
# source("../life_history_diversity_analysis/AmP_functions.R")

### AMP FUNCTIONS
getPar <- function(species, par, deb_data) {
  splist <- unlist(labels(deb_data$allStat)) #get list of all species
  par.names <- unlist(labels(deb_data$allStat[[ #extract labels for target species
    which(splist == species)
  ]]))
  par <- unlist(deb_data$allStat[[
    which(splist == species)
  ]][
    which(par.names == par)]
  )
  return(par)
}

#extracts primary climate code for given species
getPrimeClim <- function(x){
  full <- x[[1]]
  return(substring(full,1,1))
}

#extarcts full set of climate codes for given species
getClim <- function(x){
  full <- x[[1]]
  return(full)
}

getMig <- function(species, deb_data) {
  splist <- unlist(labels(deb_data$allStat)) #get list of all species
  par.names <- unlist(labels(deb_data$allStat[[ #extract labels for target species
    which(splist == species)
  ]]))
  eco <- deb_data$allStat[[
    which(splist == species)
  ]][which(par.names == "ecoCode")]
  mig <- ifelse(length(eco[[1]][[5]]) == 0, "No", eco[[1]][[5]])
  return(mig)
}

getHab <- function(species, deb_data) {
  splist <- unlist(labels(deb_data$allStat)) #get list of all species
  par.names <- unlist(labels(deb_data$allStat[[ #extract labels for target species
    which(splist == species)
  ]]))
  eco <- deb_data$allStat[[
    which(splist == species)
  ]][which(par.names == "ecoCode")]
  hab <- ifelse(length(eco[[1]][[3]]) == 0, NA, eco[[1]][[3]])
  return(hab)
}
###


deb <- allStat

#pull species list over which to iterate accessor functions
splist <- unlist(labels(deb$allStat))
splist <- splist[1:(length(splist)-2)] #remove last two entries which are not species

#extract ecoCode vector
eco <- lapply(splist, getPar, par = "ecoCode", deb_data = deb)


#pull species list over which to iterate accessor functions
splist <- unlist(labels(deb$allStat))

splist <- splist[1:(length(splist)-2)] #remove last two entries which are not species

deb <- allStat

deb_df <- data.frame(
  species = unlist(splist),
  model = unlist(lapply(splist, getPar, par ="model", deb_data = deb)), #which version of DEB model
  complete = unlist(lapply(splist, getPar, par ="COMPLETE", deb_data = deb)), #model completeness score
  climate = unlist(lapply(eco, getPrimeClim)), #primary climate code
  climate_full = unlist(lapply(eco, getClim)), #full climate code
  class = unlist(lapply(splist, getPar, par ="class", deb_data = deb)), #taxonomic class
  family = unlist(lapply(splist, getPar, par ="family", deb_data = deb)), #taxonomic family
  phylum = unlist(lapply(splist, getPar, par ="phylum", deb_data = deb)), #taxonomic phylum
  #mig = unlist(lapply(eco, getMig)), #migratory status
  mig = unlist(lapply(splist, getMig, deb_data = deb)),
  habitat = unlist(lapply(splist, getHab, deb_data = deb)),
  max_l = unlist(lapply(splist, getPar, par ="L.m", deb_data = deb)), #max length
  E_0 = unlist(lapply(splist, getPar, par ="E.0", deb_data = deb)), #embryo cost
  z = unlist(lapply(splist, getPar, par ="z", deb_data = deb)), #zoom factor
  kappa = unlist(lapply(splist, getPar, par = "kap", deb_data = deb)), #somatic allocation fraction (kappa)
  k.m = unlist(lapply(splist, getPar, par ="k.M", deb_data = deb)), #maintenance ratio (k)
  g = unlist(lapply(splist, getPar, par ="g", deb_data = deb)), #energy investment ratio (g)
  v=unlist(lapply(splist, getPar, par ="v", deb_data = deb)),
  E.m=unlist(lapply(splist, getPar, par ="E.m", deb_data = deb)),
  p.m = unlist(lapply(splist, getPar, par ="p.M", deb_data = deb)),
  l.i = unlist(lapply(splist, getPar, par ="L.i", deb_data = deb)),
  maturity_ratio = unlist(lapply(splist, getPar, par ="s.Hbp", deb_data = deb))
  ) %>% 
  mutate(
    pol = (kappa*v)/z
  )

write_csv(deb_df, "./data/deb_db_05262023.csv")


