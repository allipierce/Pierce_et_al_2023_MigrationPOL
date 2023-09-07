library(tidyverse)

##change directory as needed
datalist <- list.files("./gaemm_out", full.names = T)
datalist <- str_subset(datalist, "GAEMMrun")

data <- data.frame()
for(i in 1:length(datalist)){
  temp <- readRDS(datalist[i])
  temp$iteration <- i
  data <- bind_rows(data, temp)
}

data <-
  data %>%
  filter(alive == TRUE) %>%
  filter(fitness > 0) %>%
  group_by(migtype) %>%
  filter(fitness >= quantile(fitness, .99)) %>%
  ungroup()

#saveRDS(data, file="./gaemm_out/GAEMMdata.Rds")
#data <- readRDS(file = "./gaemm_out/GAEMMdata.Rds")

