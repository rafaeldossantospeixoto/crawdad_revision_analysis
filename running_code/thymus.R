library(tidyverse)
library(crawdad)

data_folder <- '/Users/rafaeldossantospeixoto/Library/CloudStorage/OneDrive-JohnsHopkins/jefworks/crawdad/repos/spatial_neighborhood_analysis/data'

## HBM757.VHCK.858 ## 11
t11 <- read.csv(paste0(data_folder, '/CODEX/HBM757.VHCK.858.meta.csv.gz'), row.names=1)
t11 <- t11 %>% rename(celltypes = 'com_nn50_VolnormExpr_data_annots_byXSQZ')
head(t11)

paste(min(t11$x), max(t11$x), min(t11$y), max(t11$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/053544cd63125fc25f6a71a8f444bafc
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

t11$x <- t11$x * x_res
t11$y <- t11$y * y_res
paste(min(t11$x), max(t11$x), min(t11$y), max(t11$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"



## HBM654.KTJK.968 ## 11
t21 <- read.csv(paste0(data_folder, '/CODEX/HBM654.KTJK.968.meta.csv.gz'), row.names=1)
t21 <- t21 %>% rename(celltypes = 'com_nn50_VolnormExpr_data_annots_byDMXZ')
head(t21)

paste(min(t21$x), max(t21$x), min(t21$y), max(t21$y))
## [1] "0 12095 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/827794b6592c5370801b97de0d9e2e4b
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

t21$x <- t21$x * x_res
t21$y <- t21$y * y_res
paste(min(t21$x), max(t21$x), min(t21$y), max(t21$y))
## [1] "0 4564.699519789 -3423.4302888802 0"

# Plot cts ----------------------------------------------------------------

t11 %>% 
  filter(celltypes != 'indistinct') %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes), size = .025) +
  scale_colour_manual(values = rainbow(8)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_minimal() +
  theme(legend.position="none")

t21 %>% 
  filter(celltypes != 'indistinct') %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes), size = .025) +
  scale_colour_manual(values = rainbow(8)) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_minimal() +
  theme(legend.position="none")


# CRAWDAD -----------------------------------------------------------------

cells <- crawdad::toSF(pos = t11[,c("x", "y")], celltypes = t11$celltypes)
scales <- seq(100, 3000, by=100)

shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 5,
                                            ncores = 7,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 10.02 mins

results <- crawdad::findTrends(cells,
                               dist = 50,
                               shuffle.list = shuffle.list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat11 <- crawdad::meltResultsList(results, withPerms = TRUE)
## Time was 113.36 mins
saveRDS(dat11, "running_code/processed_data/dat11.RDS")

ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

vizColocDotplot(dat11, zsig.thresh = zsig, zscore.limit = 2*zsig,
                dot.sizes = c(2,17)) +
  theme(axis.text.x = element_text(angle = 35, h = 0))



## convert dataframe to spatial points (SP)
cells <- crawdad::toSP(pos = t21[,c("x", "y")], celltypes = t21$celltypes)
## define the scales to analyze the data
scales <- c(100, 250, 500, 750, 1000, 1500)
## shuffle cells to create null background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = 7,
                                            seed = 1,
                                            verbose = TRUE)
## calculate the zscore for the cell-type pairs at different scales
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)
## summary visualization
p <- vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
p
pdf('rafael/thymus/coloc21.pdf', width = 7, height = 5)
p
dev.off()
