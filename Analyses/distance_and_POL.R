# Empirically characterizing the relationship between seasonal niche dissimilarity and pace-of-life (POL).
# SW Yanco 2023

####---- Initialization ----####

library(tidyverse)
library(ggplot2)
library(ape)
library(brms)
library(loo)
library(viridisLite)
library(ggthemes)
library(bayestestR)
library(parameters)
library(emmeans)
library(tidybayes)
library(marginaleffects)
library(patchwork)
library(performance)
library(glue)

####---- Load Data ----####

# Niche dissimilarity from Cohen et al 2023
niche0 <- read_delim("./data/niche_dissim_cohen.csv")

# Migration Distance (calculated in 'calc-dist.r')
# you can change this location to wherever you save the output from 'cal_dist.r'
dist0 <- read.csv("./out/mig_dist.csv") %>%
  distinct() # remove duplicates

# Load taxonomies
vertnet_tax <- read_csv("./data/vertlife_taxonomies.csv")

# Load AmP Data
amp0 <- read.csv("./data/deb_db_05262023.csv")

####---- Process data ----####

#Data cleaning, filtering, transforming, etc.

#Fix distance column names and add derived distance in y direction only
colnames(dist0) <- c("species", "x_breed", "y_breed", "x_wint", "y_wint", "dist")

# Add y distance
dist <- dist0 %>%
  mutate(y_dist = y_breed - y_wint) # add "y only" distance calc

# Fix species names and filter implausible POL estimates
amp <- amp0 %>%
  mutate(sp2 = str_replace(species, "\\.", " ")) %>% # fix species names
  filter(pol <1) # remove erroneous POL measures

#Combine datasets into dfs for modeling/plotting
niche_comb <- amp %>%
  inner_join(niche0, by = c("sp2" = "Species")) %>%
  inner_join(dist, by = c("sp2" = "species")) %>%
  inner_join(vertnet_tax, by = c("sp2" = "scientificname")) %>%
  mutate(sp3 = str_replace(species, "\\.", "_"),
         dist_scale = scale(dist),
         log_dist = scale(log(dist + 0.00000000001)),
         niche_scale = scale(temp_MD),
         breed_scale = scale(y_breed)) %>%
  select(pol,kappa, niche_scale, v, dist_scale,
         breed_scale,
         maturity_ratio, log_dist, sp3) %>%
  filter(complete.cases(.))


# Add phylogeny
birds <- read.nexus("./data/tree-pruner-6598654f-18d0-4125-9e50-004bf82af8de/output.nex")
birds <- birds[[runif(1, min = 1, max = length(birds))]]
phylo_vcov <- ape::vcv.phylo(birds, corr = T)

mod_dat <- niche_comb %>%
  filter(sp3 %in% rownames(phylo_vcov)) #not great but restrict to names contained in


####---- Fit Models ----####

##-- Kappa Model --##
form_kappa <- bf(pol ~ 1 + breed_scale + dist_scale + niche_scale + (1|gr(sp3, cov = phylo_vcov)),
                 phi ~ niche_scale)

mod_kappa <- brm(form_kappa,
                 # data = niche_comb,
                 data = mod_dat,
                 data2 = list(phylo_vcov = phylo_vcov),
                 family = Beta(),
                 init = 0,
                 cores = 4,
                 warmup = 5000,
                 iter = 10000,
                 thin = 5,
                 control = list(adapt_delta = .99),
                 save_pars = save_pars(all = TRUE))


# Model check
pp_check(mod_kappa)


##-- Maturity Model --##
form_mat <- bf(maturity_ratio ~ 1 + breed_scale + dist_scale + niche_scale + (1|gr(sp3, cov = phylo_vcov)))

mod_mat <- brm(form_mat,
               # data = niche_comb,
               data = mod_dat,
               data2 = list(phylo_vcov = phylo_vcov),
               family = Beta(),
               init = 0,
               cores = 4,
               warmup = 5000,
               iter = 50000,
               thin = 7,
               control = list(adapt_delta = .99),
               save_pars = save_pars(all = TRUE))


# Check Models
pp_check(mod_mat)


## Phylogenetic Signal
(vc_kappa <- variance_decomposition(mod_kappa, robust = T))
(vc_mat <- variance_decomposition(mod_mat,robust = T))


# Save Models
save(mod_kappa, mod_mat, file = glue("../out/models_{Sys.Date()}.rdata"))

# Load Models
load("./out/models_2023-06-13.rdata")

####---- Plot models ----####

##-- Kappa Component Plots --##

# NICHE
trends_kappa_niche <- mod_kappa %>%
  emtrends(~ niche_scale,
           var = "niche_scale",
           at = list(niche_scale = c(0)),
           type = "response") %>%
  gather_emmeans_draws() %>%
  mutate(par = "niche")

(pd_kappa_niche <- sum(trends_kappa_niche$.value > 0)/nrow(trends_kappa_niche))
median(trends_kappa_niche$.value)

(kappa_niche_plot <- ggplot(trends_kappa_niche, aes(x = .value, fill = factor(niche_scale))) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.95), point_interval = "median_hdi",
                 slab_alpha = 0.75) +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    labs(x = "Niche Dissimilarity",
         y = "Density", fill = "niche_scale",
         caption = glue("pd = {round(pd_kappa_niche*100, 2)}%")) +
    theme_clean() +
    theme(legend.position = "none",
          text=element_text(family="Arial")))

# DISTANCE
trends_kappa_dist <- mod_kappa %>%
  emtrends(~ dist_scale, var = "dist_scale",
           at = list(dist_scale = c(0)),
           regrid = "response") %>%
  gather_emmeans_draws() %>%
  mutate(par = "dist")

(pd_kappa_dist <- sum(trends_kappa_dist$.value > 0)/nrow(trends_kappa_dist))
median(trends_kappa_dist$.value)

(kappa_dist_plot <- ggplot(trends_kappa_dist, aes(x = .value, fill = factor(dist_scale))) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.95), point_interval = "median_hdi",
                 slab_alpha = 0.75) +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    labs(x = "Migration Distance",
         y = "Density", fill = "niche_scale",
         caption = glue("pd = {round(pd_kappa_dist*100, 2)}%")) +
    theme_clean() +
    theme(legend.position = "none",
          text=element_text(family="Arial")))

# BREEDING LOCATION
trends_kappa_breed <- mod_kappa %>%
  emtrends(~ breed_scale, var = "breed_scale",
           at = list(breed_scale = c(0)),
           regrid = "response") %>%
  gather_emmeans_draws() %>%
  mutate(par = "breed")

(pd_kappa_breed <- sum(trends_kappa_breed$.value > 0)/nrow(trends_kappa_breed))
median(trends_kappa_breed$.value)

(kappa_breed_plot <- ggplot(trends_kappa_breed, aes(x = .value, fill = factor(breed_scale))) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.95), point_interval = "median_hdi",
                 slab_alpha = 0.75) +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    labs(x = "Breeding Latitude",
         y = "Density", fill = "niche_scale",
         caption = glue("pd = {round(pd_kappa_breed*100, 2)}%")) +
    theme_clean() +
    theme(legend.position = "none",
          text=element_text(family="Arial")))

##-- Maturity Ratio Component Plots --##

#NICHE
trends_mat_niche <- mod_mat %>%
  emtrends(~ niche_scale,
           var = "niche_scale",
           at = list(niche_scale = c(0)),
           type = "response") %>%
  gather_emmeans_draws()

(pd_mat_niche <- sum(trends_mat_niche$.value > 0)/nrow(trends_mat_niche))
median(trends_mat_niche$.value)

(mat_niche_plot <- ggplot(trends_mat_niche, aes(x = .value, fill = factor(niche_scale))) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.95), point_interval = "median_hdi",
                 slab_alpha = 0.75) +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    labs(x = "Niche Dissimilarity",
         y = "Density", fill = "niche_scale",
         caption = glue("pd = {round(pd_mat_niche*100, 2)}%")) +
    theme_clean() +
    theme(legend.position = "none",
          text=element_text(family="Arial")))

# DISTANCE
trends_mat_dist <- mod_mat%>%
  emtrends(~ dist_scale, var = "dist_scale",
           at = list(dist_scale = c(0)),
           regrid = "response") %>%
  gather_emmeans_draws()

(pd_mat_dist <- sum(trends_mat_dist$.value > 0)/nrow(trends_mat_dist))
median(trends_mat_dist$.value)

(mat_dist_plot <- ggplot(trends_mat_dist, aes(x = .value, fill = factor(dist_scale))) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.95), point_interval = "median_hdi",
                 slab_alpha = 0.75) +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    labs(x = "Migration Distance",
         y = "Density", fill = "niche_scale",
         caption = glue("pd = {round(pd_mat_dist*100, 2)}%")) +
    theme_clean() +
    theme(legend.position = "none",
          text=element_text(family="Arial")))

#BREEDINMG LOCATION
trends_mat_breed <- mod_mat %>%
  emtrends(~ breed_scale, var = "breed_scale",
           at = list(breed_scale = c(0)),
           regrid = "response") %>%
  gather_emmeans_draws()

(pd_mat_breed <- sum(trends_mat_breed$.value < 0)/nrow(trends_mat_breed))
median(trends_mat_breed$.value)

(mat_breed_plot <- ggplot(trends_mat_breed, aes(x = .value, fill = factor(breed_scale))) +
    geom_vline(xintercept = 0) +
    stat_halfeye(.width = c(0.95), point_interval = "median_hdi",
                 slab_alpha = 0.75) +
    scale_fill_viridis_d(option = "viridis", end = 0.6) +
    labs(x = "Breeding Latitude",
         y = "Density", fill = "niche_scale",
         caption = glue("pd = {round(pd_mat_breed*100, 2)}%")) +
    theme_clean() +
    theme(legend.position = "none",
          text=element_text(family="Arial")))


##-- Combine plots --##

(comb_plot <- wrap_elements(kappa_niche_plot/kappa_dist_plot/kappa_breed_plot+plot_annotation(title = "POL"))|wrap_elements(mat_niche_plot/mat_dist_plot/mat_breed_plot+plot_annotation(title = "Maturity Ratio")))

ggsave(comb_plot, filename = "../out/empirical.png", dpi = 300, width = 8, height = 6)
ggsave(plot = comb_plot, filename = "../out/empirical.pdf", width = 8, height = 6, device = cairo_pdf)

