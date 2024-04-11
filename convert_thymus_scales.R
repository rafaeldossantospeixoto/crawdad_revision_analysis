library(crawdad)
library(tidyverse)

# convert vhck ------------------------------------------------------------

## HBM757.VHCK.858 ## 11 yrs
vhck <- read.csv('../spatial_neighborhood_analysis/data/CODEX/HBM757.VHCK.858.meta.csv.gz', row.names=1)
vhck <- vhck %>% 
  rename(celltypes = 'com_nn50_VolnormExpr_data_annots_byXSQZ') %>% 
  select(c('x', 'y', 'celltypes')) %>% 
  dplyr::mutate_at(vars(x, y), as.numeric)
head(vhck)

paste(min(vhck$x), max(vhck$x), min(vhck$y), max(vhck$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/053544cd63125fc25f6a71a8f444bafc
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

vhck$x <- vhck$x * x_res
vhck$y <- vhck$y * y_res
paste(min(vhck$x), max(vhck$x), min(vhck$y), max(vhck$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

write.csv2(vhck, file = paste0('running_code/data/thymus/VHCK.meta.csv.gz'), row.names = TRUE)



# convert ktjk ------------------------------------------------------------

## HBM654.KTJK.968 ## 21 yrs
ktjk <- read.csv('../spatial_neighborhood_analysis/data/CODEX/HBM654.KTJK.968.meta.csv.gz', row.names=1)
ktjk <- ktjk %>% 
  rename(celltypes = 'com_nn50_VolnormExpr_data_annots_byDMXZ') %>% 
  select(c('x', 'y', 'celltypes')) %>% 
  dplyr::mutate_at(vars(x, y), as.numeric)
head(ktjk)

paste(min(ktjk$x), max(ktjk$x), min(ktjk$y), max(ktjk$y))
## [1] "0 12095 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/827794b6592c5370801b97de0d9e2e4b
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

ktjk$x <- ktjk$x * x_res
ktjk$y <- ktjk$y * y_res
paste(min(ktjk$x), max(ktjk$x), min(ktjk$y), max(ktjk$y))
## [1] "0 4564.699519789 -3423.4302888802 0"

write.csv2(ktjk, file = paste0('running_code/data/thymus/KTJK.meta.csv.gz'), row.names = TRUE)

