library(data.table)
library(doRNG)
library(furrr)
library(tidyverse)

source("R/expendenergy.R")
source("R/move.R")
source("R/reproduce.R")
source("R/randomizetraits.R")
source("R/GAEMM.R")
source("R/seasonalworld.R")
source("R/egg_vals.R")
source("R/checkstates.R")


basemodelparams <-  list(worldsize = 180,
                         popsize = 5000,
                         steps = 365,
                         generations = 200,
                         kappa_min = 0.01,
                         kappa_max = 0.99,
                         maxmoves = 45,
                         e_traits = c("kappa", "m", "v", "g", "z", "l_b","uph"),
                         m_traits = c("mig_time", "nb_length", "a_b"),
                         startloc = 80,
                         endloc = 80,
                         replacerate = 1,
                         crossoverrate = 0.5,
                         mutationrate = 0.5,
                         selectiontype = "sus",
                         niching = T,
                         movement = F,
                         meanrange = c(0.6,0.5))
save(basemodelparams,file="tuningdash/basemodelparams.rda")
reps <- 5

## tuning param conditions
popsizes <- c(2500, 5000, 10000)
generations <- 500
replacerates <- c(0.5, 0.75, 1)
crossoverrates <- c(0.25, 0.5, 0.75)
mutationrates <- c(0.1, 0.25, 0.5)

modparamlist <- function(x, paramname, modelparams) {
  modelparams
  modelparams[[paramname]] <- x
  modelparams$paramname <- paramname
  modelparams
}

ps <- lapply(popsizes, FUN = modparamlist, paramname = "popsize", modelparams = basemodelparams)
gt <- lapply(generations, FUN = modparamlist, paramname = "generation", modelparams = basemodelparams)
rr <- lapply(replacerates, FUN = modparamlist, paramname = "replacerate", modelparams = basemodelparams)
cr <- lapply(crossoverrates, FUN = modparamlist, paramname = "crossoverrate", modelparams = basemodelparams)
mr <- lapply(mutationrates, FUN = modparamlist, paramname = "mutationrate", modelparams = basemodelparams)

modelparams <- c(ps,rr,cr,mr,gt)
modelparams <- rep(modelparams, reps)
i <- rep(1:reps, each = length(modelparams))
modelparams <- Map(c, modelparams, i = i)
seed <- 42

## This is set up so that random number generation is independent between iterations but is fixed across movement conditions, i.e. movement conditions are run with the same starting populations (and new random numbers are generated from the same seeds).

runtuneGAEMM <- function(x){
  data.table::setDTthreads(1)
  data <- GAEMM(modelparams = x, lastonly = T, verbose = F, tuning = T)
  data$paramname <- x$paramname
  data$paramvalue <- x[[x$paramname]]
  data$ids <- 1:nrow(data)
  saveRDS(data, file=paste0("tuningdat_",x$paramname,"_",x[[x$paramname]],"_",x$i,".Rds"))
  rm(data)
  gc()
}

plan(multisession, workers = 10, gc = TRUE)

rng <- RNGseq(reps, seed)
rng <- rep(rng, each = length(modelparams)/reps)

#future.scheduling false makes each future start a fresh instance of r
#future_lapply(modelparams, runparGAEMM, future.seed = rng, future.packages = "data.table", future.scheduling = FALSE)

future_walk(modelparams, runtuneGAEMM, .options = furrr_options(seed = rng, stdout = F))
plan(sequential)

datalist <- list.files(".", full.names = T)
datalist <- str_subset(datalist, "tuningdat")
tuningdat <-  datalist %>%
  map_dfr(readRDS)
save(tuningdat,file="tuningdash/tuningdat.rda")
file.remove(datalist)
rm(tuningdat)
rm(datalist)
