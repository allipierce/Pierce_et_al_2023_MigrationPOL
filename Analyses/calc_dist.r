# Calculate range centroids and distance for empirical data analysis
# SW Yanco 2023

library(sf)
library(tidyverse)
# library(geosphere)

dat0 <- st_read("data/BOTW/BOTW.gdb")
file.create("out/mig_dist.csv")

sp_l <- dat0 %>%
  select(SCINAME) %>%
  unique()
st_geometry(sp_l) <- NULL

for (i in 1:nrow(sp_l)){
  dat1 <- dat0 %>%
    filter(SCINAME == sp_l[i,])

  sp <- sp_l[i,]

  b <- st_combine(bind_rows(dat1[dat1$SEASONAL == 1,], dat1[dat1$SEASONAL == 2,])) %>%
    st_transform(5070) %>%
    st_centroid()

  w <- st_combine(bind_rows(dat1[dat1$SEASONAL == 3,], dat1[dat1$SEASONAL == 1,])) %>%
    st_transform(5070) %>%
    st_centroid()

  d <- st_distance(b, w)

  out <- matrix(data = c(sp,
                         st_coordinates(b)[1],
                         st_coordinates(b)[2],
                         st_coordinates(w)[1],
                         st_coordinates(w)[2],
                         d),
                nrow = 1)
  write.table(out, file = "out/mig_dist.csv", append = T,row.names = F, col.names = F, sep = ",")
}
