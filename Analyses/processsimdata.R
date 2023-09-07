## Simulation data processing/loading
# follows directory structure output of parallelized model scripts

library(tidyverse)
library(data.table)
source("./R/seasonalworld.R")

# load minimal model params for data processing
modelparams <- list(worldsize = 180,
                    popsize = 2000,
                    steps = 365,
                    e_traits = c("kappa", "m", "v", "g", "z", "l_b","uph"),
                    m_traits = c("mig_time", "nb_length", "a_b"),
                    meanrange = c(0.6,0.5))

utraits <- c("kappa", "m", "g", "l_b", "l_p", "v", "mig_time", "nb_length", "a_b", "fitness")

### filter data from 1000 runs to top 5 results per movement strategy

 path <- "./data"
 dirlist <- list.dirs(path, full.names = T)
 dirlist <- str_subset(dirlist, "seed")
 datalist <-
   unlist(lapply(
     dirlist,
     list.files,
     full.names = T,
     pattern = "GAEMMrun"
   ))

 migtypes <- c("LD_HS","LD_MS","MD_AS","MD_MS","SD_MS","SD_AS","NM_HS","NM_MS","NM_AS")
 names(migtypes) <- migtypes

 datalists <-
   lapply(migtypes, function(x)
     str_subset(datalist, pattern = x))

top5result <- function(f) {
  unique(as.data.table(read_rds(f)), by = c(utraits)) %>%
    filter(alive == TRUE) %>%
    filter(fitness > 0) %>%
    mutate(across(c(kappa,m,v,g,z,l_b,l_p,mig_time,nb_length,a_b), round, 5)) %>%
    distinct(kappa,m,v,g,z,l_b,l_b,mig_time,nb_length,a_b,migtype, .keep_all = T) %>%
    slice_max(fitness, n = 5, with_ties = T)
}

 initdata <- as.list(rep(NA, length(migtypes)))
 names(initdata) <- migtypes

 for (i in 1:length(migtypes)) {
   initdata[[i]] <- map_dfr(datalists[[i]],top5result,.progress = T)
 }

 saveRDS(initdata, file = "./data/initdata.Rds")

### load and filter data from pre-initialized runs

datalist <-
  list.files("./data/rerunpreinit/",
             full.names = T,
             pattern = "GAEMMtoprun")
 datalists <-
   lapply(migtypes, function(x)
     str_subset(datalist, pattern = x))

rerundata <- as.list(rep(NA, length(migtypes)))
for (i in 1:length(migtypes)) {
  rerundata[[i]] <- map_dfr(datalists[[i]],top5result,.progress = T)
}
rerundata <- bind_rows(rerundata)

rerundata <- rerundata %>%
  filter(alive == TRUE) %>%
  filter(fitness > 0) %>%
  distinct(kappa,
           m,
           v,
           g,
           z,
           l_b,
           l_p,
           mig_time,
           nb_length,
           a_b,
           migtype,
           .keep_all = T)

rerundata_fil <-
  rerundata %>%
  group_by(migtype) %>%
  slice_max(fitness, n = 5) %>%
  ungroup()

rerundata_fil <-
  rerundata_fil %>%
  mutate(mig_time = case_match(migtype,
                               c("NM_AS", "NM_MS", "NM_HS") ~ 0,
                               .default = mig_time)) %>%
  mutate(nb_length = case_match(migtype,
                                c("NM_AS", "NM_MS", "NM_HS") ~ 0,
                                .default = nb_length)) %>%
  mutate(m = case_match(migtype,
                        c("NM_AS", "NM_MS", "NM_HS") ~ 0,
                        .default = m)) %>%
  mutate(mig_speed = case_match(migtype,
                                c("NM_AS", "NM_MS", "NM_HS") ~ 0,
                                .default = mig_speed))


## Combine data, filter, and save
data <- bind_rows(read_rds("./data/initdata.Rds"))
data <- bind_rows(data, rerundata_fil)

# filter and annotate combined data for plotting/analysis
df <- data %>%
  filter(alive == TRUE) %>%
  filter(fitness > 0) %>%
  distinct(kappa,
           m,
           v,
           g,
           z,
           l_p,
           l_b,
           mig_time,
           nb_length,
           a_b,
           migtype,
           .keep_all = T)

df <- df %>%
  mutate(
    movetype = case_when(
      migtype == "LD_HS" ~ "M",
      migtype == "LD_MS" ~ "M",
      migtype == "MD_AS" ~ "M",
      migtype == "MD_MS" ~ "M",
      migtype == "SD_AS" ~ "M",
      migtype == "SD_MS" ~ "M",
      migtype == "NM_HS" ~ "NM",
      migtype == "NM_MS" ~ "NM",
      migtype == "NM_AS" ~ "NM"
    )
  ) %>%
  mutate(bseason = case_when(startloc == 80 ~ "HS",
                             startloc == 40 ~ "MS",
                             startloc == 0 ~ "AS")) %>%
  mutate(
    wseason = case_when(
      endloc == -80 ~ "HS",
      endloc == 80 ~ "HS",
      endloc == -40 ~ "MS",
      endloc == 40 ~ "MS",
      endloc == 0 ~ "AS"
    )
  ) %>%
  mutate(
    mdist = case_when(
      migtype == "LD_HS" ~ "LD",
      migtype == "LD_MS" ~ "LD",
      migtype == "MD_AS" ~ "MD",
      migtype == "MD_MS" ~ "MD",
      migtype == "SD_AS" ~ "SD",
      migtype == "SD_MS" ~ "SD",
      migtype == "NM_HS" ~ "NM",
      migtype == "NM_MS" ~ "NM",
      migtype == "NM_AS" ~ "NM"
    )
  ) %>%
  mutate(
    subtype = case_when(
      migtype == "LD_HS" ~ "a",
      migtype == "LD_MS" ~ "b",
      migtype == "MD_AS" ~ "a",
      migtype == "MD_MS" ~ "b",
      migtype == "SD_AS" ~ "b",
      migtype == "SD_MS" ~ "a",
      migtype == "NM_HS" ~ "a",
      migtype == "NM_MS" ~ "a",
      migtype == "NM_AS" ~ "a"
    )
  ) %>%
  mutate(
    mdist = case_when(
      migtype == "LD_HS" ~ "LD",
      migtype == "LD_MS" ~ "LD",
      migtype == "MD_AS" ~ "MD",
      migtype == "MD_MS" ~ "MD",
      migtype == "SD_AS" ~ "SD",
      migtype == "SD_MS" ~ "SD",
      migtype == "NM_HS" ~ "NM",
      migtype == "NM_MS" ~ "NM",
      migtype == "NM_AS" ~ "NM"
    )
  )  %>%
  mutate(distsub = paste(mdist, subtype)) %>%
  mutate(fdiff = abs(
    lookupf_bylat(1, startloc, modelparams) -
      lookupf_bylat(181, endloc, modelparams)
  )) %>%
  mutate_if(is_character, as_factor) %>%
  mutate(migtype = fct_reorder2(migtype, endloc, mdist)) %>%
  mutate(dist = abs(startloc - endloc)) %>%
  mutate(etype = case_match(
    migtype,
    c("LD_HS", "LD_MS", "MD_MS") ~ "Buffer",
    c("MD_AS", "SD_MS", "SD_AS") ~ "Mitigate",
    c("NM_AS", "NM_MS", "NM_HS") ~ "Cope"
  ))

saveRDS(df, "./data/df.Rds")

mergesimilar <-
  df %>%
  mutate(ID = 1:nrow(df)) %>%
  mutate(across(all_of(
    c(
      "kappa",
      "m",
      "v",
      "g",
      "z",
      "l_p",
      "l_b",
      "mig_time",
      "nb_length",
      "a_b"
    )
  ), round, digits = 5)) %>%
  distinct(kappa,
           m,
           v,
           g,
           z,
           l_p,
           l_b,
           mig_time,
           nb_length,
           a_b,
           migtype,
           .keep_all = T) %>%
  pull(ID)

dftop <-
  df[mergesimilar, ] %>%
  mutate(
    mig_speed = case_when(
      migtype == "NM_HS" ~ 0,
      migtype == "NM_MS" ~ 0,
      migtype == "NM_AS" ~ 0,
      TRUE ~ mig_speed,
    )
  ) %>%

  group_by(migtype) %>%
  filter(fitness / max(fitness) >= .99) %>%
  ungroup()

dftop$ID <- 1:nrow(dftop)
saveRDS(dftop, "./data/dftop.Rds")

# create to long form version data for path plotting
dfsteps <- dftop %>%
  dplyr::select(!loc) %>%
  pivot_longer(
    cols = starts_with("step"),
    names_to = "step",
    names_prefix = "step",
    values_to = "loc"
  ) %>%
  mutate(step = as.integer(step))

dfsteps <- dfsteps %>%
  mutate(f = lookupf_bylat(step, loc, modelparams)) %>%
  mutate(ss = lookupf_bylat(step, startloc, modelparams)) %>%
  mutate(es = lookupf_bylat(step, endloc, modelparams))

saveRDS(dfsteps, "./data/dfsteps.Rds")
