# Script to run pre-initialized parameters in parallel for each movement strategy

library(data.table)
library(tidyverse)
library(doRNG)
library(furrr)

source("R/expendenergy.R")
source("R/move.R")
source("R/reproduce.R")
source("R/randomizetraits.R")
source("R/GAEMM.R")
source("R/seasonalworld.R")
source("R/egg_vals.R")
source("R/checkstates.R")

modelparams <- list()
modelparams[[1]] <- list(worldsize = 180,
                         popsize = 5000,
                         steps = 365,
                         generations = 200,
                         kappa_min = 0.01,
                         kappa_max = 0.99,
                         maxmoves = 45,
                         e_traits = c("kappa", "m", "v", "g", "z", "l_b","uph"),
                         m_traits = c("mig_time", "nb_length", "a_b"),
                         startloc = 80,
                         endloc = -80,
                         replacerate = 1,
                         crossoverrate = 0.5,
                         mutationrate = 0.5,
                         selectiontype = "sus",
                         niching = T,
                         movement = T,
                         meanrange = c(0.6,0.5))

names(modelparams)[1] <- "LD_HS"
modelparams[[1]]$migtype <- "LD_HS"

modelparams[[2]] <- modelparams[[1]]
names(modelparams)[2] <- "LD_MS"
modelparams[[2]]$endloc <- -40
modelparams[[2]]$migtype <- "LD_MS"

modelparams[[3]] <- modelparams[[2]]
names(modelparams)[3] <- "MD_AS"
modelparams[[3]]$endloc <- 0
modelparams[[3]]$migtype <- "MD_AS"

modelparams[[4]] <- modelparams[[3]]
names(modelparams)[4] <- "MD_MS"
modelparams[[4]]$startloc <- 40
modelparams[[4]]$endloc <- -40
modelparams[[4]]$migtype <- "MD_MS"

modelparams[[5]] <- modelparams[[4]]
names(modelparams)[5] <- "SD_MS"
modelparams[[5]]$startloc <- 80
modelparams[[5]]$endloc <- 40
modelparams[[5]]$migtype <- "SD_MS"

modelparams[[6]] <- modelparams[[5]]
names(modelparams)[6] <- "SD_AS"
modelparams[[6]]$startloc <- 40
modelparams[[6]]$endloc <- 0
modelparams[[6]]$migtype <- "SD_AS"

modelparams[[7]] <- modelparams[[1]]
names(modelparams)[7] <- "NM_HS"
modelparams[[7]]$movement <- FALSE
modelparams[[7]]$startloc <- 80
modelparams[[7]]$endloc <- 80
modelparams[[7]]$migtype <- "NM_HS"

modelparams[[8]] <- modelparams[[7]]
names(modelparams)[8] <- "NM_MS"
modelparams[[8]]$startloc <- 40
modelparams[[8]]$endloc <- 40
modelparams[[8]]$migtype <- "NM_MS"

modelparams[[9]] <- modelparams[[8]]
names(modelparams)[9] <- "NM_AS"
modelparams[[9]]$startloc <- 0
modelparams[[9]]$endloc <- 0
modelparams[[9]]$migtype <- "NM_AS"

utraits <- c("kappa", "m", "g", "l_b", "l_p", "v", "mig_time", "nb_length", "a_b", "fitness")

### Prep initial data for first iteration
 # initdata <- readRDS(file = "./data/initdata.Rds")
 # dir.create(file.path("./data/initdata"), showWarnings = FALSE)
 # walk2(initdata, names(initdata), function(x, y) saveRDS(x, file = paste0("./data/initdata/initdata_",y,"_001.Rds")))

# Function to filter results before writing to disk
top5result_df <- function(df){
  unique(as.data.table(df), by = c(utraits)) %>%
    slice_max(fitness, n = 5, with_ties = T)
}

reps <- 10
seed <- 6

runparGAEMM <- function(x,y){
  data.table::setDTthreads(1)
  idata <- read_rds(y)
  data <- GAEMM(modelparams = x, lastonly = T, verbose = F, initdata = idata)
  data$migtype <- x$migtype
  data$i <- x$i
  data$seedreps <- paste0(seed,"_",reps)
  saveRDS(data, file=paste0("./data/rerunpreinit/GAEMMtoprun_",sprintf("%03d", x$i),"_",x$migtype,".Rds"))
  saveRDS(top5result_df(data), file = paste0("./data/initdata/initdata_",x$migtype,"_",sprintf("%03d", x$i + 1),".Rds"))
  rm(data)
  rm(idata)
  gc()
}


plan(multisession, workers = 9, gc = TRUE)

## This is set up so that random number generation is independent between iterations but is fixed across movement conditions, i.e. movement conditions are run with the same starting populations (and new random numbers are generated from the same seeds).
rng <- RNGseq(reps, seed)

for(i in 1:reps){
  plan(multisession, workers = 9, gc = TRUE)
  initdatanames <- paste0("./data/initdata/initdata_",names(modelparams),"_",sprintf("%03d", i),".Rds")
  modelparams <- lapply(modelparams, function(x) {x$i = i; return(x)})
  future_walk2(modelparams, initdatanames, runparGAEMM, .options = furrr_options(seed = rep(rng[i], length(modelparams)), stdout = F))
  plan(sequential)
}





