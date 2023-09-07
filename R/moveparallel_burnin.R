library(data.table)
#library(doFuture)
library(tidyverse)
library(doRNG)
#library(future.apply)
library(furrr)

source("R/expendenergy.R")
source("R/move.R")
source("R/reproduce.R")
source("R/randomizetraits.R")
source("R/GAEMM_new.R")
source("R/seasonalworld.R")
source("R/egg_vals.R")
source("R/checkstates.R")

modelparams <- list()
modelparams[[1]] <- list(worldsize = 180,
                    popsize = 5000,
                    steps = 365,
                    generations = 100,
                    kappa_min = 0.01,
                    kappa_max = 0.99,
                    maxmoves = 45,
                    e_traits = c("kappa", "m", "v", "g", "z", "l_b","uph"),
                    m_traits = c("mig_time", "nb_length", "a_b"),
                    startloc = 80,
                    endloc = -80,
                    replacerate = 1,
                    crossoverrate = 0.5,
                    mutationrate = 0.05,
                    selectiontype = "sus",
                    niching = T,
                    movement = T,
                    meanrange = c(0.6,0.5),
                    ebtol = 0.001)

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

reps <- 100
i <- rep(1:reps, each = length(modelparams))

seed <- 123

dir.create(file.path("data_eb"), showWarnings = FALSE)
dir.create(file.path("./data_eb", paste0("seed",seed)), showWarnings = FALSE)

modelparams <- rep(modelparams, reps)
modelparams <- Map(c, modelparams, i = i)

## This is set up so that random number generation is independent between iterations but is fixed across movement conditions, i.e. movement conditions are run with the same starting populations (and new random numbers are generated from the same seeds).

# function to run model and save output in parallel
runparGAEMM <- function(x){
  data.table::setDTthreads(1)
  data <- GAEMM(modelparams = x, lastonly = T, verbose = F)
  data$migtype <- x$migtype
  data$i <- x$i
  data$seedreps <- paste0(seed,"_",reps)
  saveRDS(data, file=paste0("./data_eb/seed",seed,"/GAEMMrun_",sprintf("%03d", x$i),"_",x$migtype,".Rds"))
  rm(data)
  gc()
}

#treats warnings like errors - used so model got interrupted if something wonky happened
#options(future.globals.onReference = "error")
#options(warn=0)
 # data <- lapply(modelparams, runparGAEMM) %>%
 #   bind_rows(.id = "migtype")

#registerDoFuture()
#gc = TRUE cleans stuff up between futures
plan(multisession, workers = 18, gc = TRUE)
#plan(future.callr::callr, workers = 9, gc = TRUE)

rng <- RNGseq(reps, seed)
rng <- rep(rng, each = length(modelparams)/reps)

#future.scheduling false makes each future start a fresh instance of r
#future_lapply(modelparams, runparGAEMM, future.seed = rng, future.packages = "data.table", future.scheduling = FALSE)

future_walk(modelparams, runparGAEMM, .options = furrr_options(seed = rng, stdout = F))
plan(sequential)

for(n in 1:reps){
  datalist <- list.files(paste0("./data_eb/seed",seed), full.names = T)
  datalist <- str_subset(datalist, paste0("GAEMMrun_",sprintf("%03d", n)))
  gdata <-  datalist %>%
    map_dfr(readRDS) %>%
    filter(alive == TRUE) %>%
    filter(fitness > 0) %>%
    mutate(across(c(kappa,m,v,g,z,l_b,l_p,mig_time,nb_length,a_b), round, 5)) %>%
    distinct(kappa,m,v,g,z,l_b,l_b,mig_time,nb_length,a_b,migtype, .keep_all = T)
    #distinct(kappa, m,v,g,z,l_b,uph,mig_time,nb_length,a_b,migtype, .keep_all = T)
  saveRDS(gdata,file=paste0("./data_eb/seed",seed,"/DataGAEMM_",format(Sys.time(), "%b%d%_%H:%M"),"_Seed",seed,"_i",n,".Rds"))
  #file.remove(datalist)
  rm(gdata)
  rm(datalist)
}




